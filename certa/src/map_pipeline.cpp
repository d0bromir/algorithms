#include "map_pipeline.h"

#include <algorithm>
#include <atomic>
#include <stdexcept>
#include <thread>

#include "cli.h"

namespace certa {

namespace {

// Calls work(w) for w in [0, n): w = 0 on this thread, the others on threads
// of their own.
template <class F>
void run_workers(int n, F&& work) {
  std::vector<std::thread> pool;
  for (int w = 1; w < n; ++w) pool.emplace_back([&work, w] { work(w); });
  work(0);
  for (auto& th : pool) th.join();
}

// Joins a thread even if mapping throws, so the error is reported instead of
// std::terminate on a joinable thread.
struct Joiner {
  std::thread& t;
  ~Joiner() {
    if (t.joinable()) t.join();
  }
};

// Which reads each pass re-runs (the same on either index).
// Passes 2 and V: every read still uncertified, except those too short or
// too long to be split into parts.
bool retry(const Result& r) { return r.certified == 0 && r.reason != kBadLength; }

// Pass 3: reads that failed because the budget limited how many parts they
// could search.
bool retry_pass3(const Result& r, const Params& p3) {
  const int want = r.parts < p3.k + p3.min_support ? r.parts : p3.k + p3.min_support;
  return r.certified == 0 && (r.reason == kRadiusNegative || r.reason == kNotFound) && r.used < want;
}

// Pass L: reads with no locus within the radius, or a cluster too wide.
bool retry_local(const Result& r) {
  return r.certified == 0 && (r.reason == kNotFound || r.reason == kClusterTooWide);
}

// The pair pass: pairs with one certified mate and an uncertified or tied
// mate. Pairs with both mates uncertified are skipped: on HG002 under 1% of
// them get a pair certificate, at a fifth of the pass's cost.
bool wants_pair_pass(const Result& x, const Result& y) {
  const bool one = x.certified == 1 || y.certified == 1;
  return one && (x.certified != 1 || y.certified != 1 || x.n_best > 1 || y.n_best > 1);
}

// Tier R (interval mode): every read not certified once the contig check
// applies, except those too short or too long.
bool wants_tier_r(const Result& r0, const Reference& ref) {
  Result r = r0;
  check_contig(r, ref);
  return r.certified != 1 && r.reason != kBadLength;
}

// The reason an interval-mode read is counted under: kOk if certified, else
// why not (a read placed by tier S2 is counted as not found).
int interval_reason(const Result& r) {
  return r.certified == 1 ? int(kOk) : (r.reason == kOk ? int(kNotFound) : int(r.reason));
}

}  // namespace

// Per-thread scratch of the fused stage.
struct MapPipeline::ChunkScratch {
  explicit ChunkScratch(bool index2) : wx(index2 ? new Workspace : nullptr) { civ.reserve(kChunk); }
  std::unique_ptr<Workspace> ws{new Workspace}, w1{new Workspace}, w2{new Workspace}, wr{new Workspace};
  std::unique_ptr<Workspace> wx;  // the second index's passes
  PairScratch sc, psc;            // the pair pass; interval mode's mate rescue
  RepeatScratch rsc;
  std::vector<RepeatOut> civ;     // the chunk's tier R outputs
  size_t idx[kChunk], id2[kChunk], sub[kChunk];
  std::unique_ptr<Result[]> rr{new Result[kChunk]}, rx{new Result[kChunk]};
};

MapPipeline::MapPipeline(const MapConfig& c, const Reference& ref, const IndexView& view, const IndexView* view2,
                         FastqInput& reads, FastqInput* mates, const MapOutputs& out)
    : c_(c), ref_(ref), view_(view), view2_(view2), reads_(reads), mates_(mates), out_(out),
      writer_(ref, view, c.out) {}

PassCounts MapPipeline::counts() const {
  PassCounts n = gpu_counts_;
  n.add(cpu_counts_);
  return n;
}

// Reads and encodes the next batch into a slot (empty at the end of input).
// An error is kept for run() to throw.
void MapPipeline::load(int slot) {
  Slot& S = slots_[slot];
  try {
    auto t = Clock::now();
    if (!c_.paired) {
      reads_.next_batch(S.recs, c_.batch_size, c_.io_threads);
    } else {
      // Mates interleaved: read 2i is mate 1, read 2i + 1 its mate 2, each
      // parsed straight into its place (strings reused from the slot's
      // previous batch).
      FastqSpans s1, s2;
      const size_t half = std::max<size_t>(1, c_.batch_size / 2);
      const size_t n1 = reads_.next_spans(s1, half, c_.io_threads);
      const size_t n2 = mates_->next_spans(s2, half, c_.io_threads);
      if (n1 != n2) throw std::runtime_error("mate files have different numbers of reads");
      std::vector<FastqRecord>& out = S.recs;
      out.resize(2 * n1);
      parse_spans(s1, out.data(), 2, c_.io_threads);
      parse_spans(s2, out.data() + 1, 2, c_.io_threads);
      for (size_t i = 0; i < n1; ++i)
        if (out[2 * i].name != out[2 * i + 1].name)
          throw std::runtime_error("mate names differ: " + out[2 * i].name + " / " + out[2 * i + 1].name);
    }
    encode_batch(S.recs, S.batch, c_.io_threads);
    time_.read_encode += secs(t);
  } catch (...) {
    S.recs.clear();
    read_error_ = std::current_exception();
  }
}

// Re-runs reads `redo` of a batch with p, as a compacted batch; with
// only_certified, a result replaces the old one only when certified. Returns
// the results taken.
size_t MapPipeline::rerun(Slot& S, const std::vector<size_t>& redo, const Params& p, bool on_gpu,
                          bool only_certified) {
  const ReadBatch& batch = S.batch;
  if (redo.empty()) return 0;
  ReadBatch sub;
  sub.offs.resize(redo.size());
  sub.lens.resize(redo.size());
  sub.hashes.resize(redo.size());
  for (size_t x = 0; x < redo.size(); ++x) {
    const size_t i = redo[x];
    sub.offs[x] = sub.codes.size();
    sub.lens[x] = batch.lens[i];
    sub.hashes[x] = batch.hashes[i];
    sub.codes.insert(sub.codes.end(), batch.codes.begin() + batch.offs[i],
                     batch.codes.begin() + batch.offs[i] + batch.lens[i]);
  }
  std::vector<Result> res2;
  if (on_gpu) gm_->map(p, sub, res2);
  else map_cpu(view_, p, sub, res2, c_.threads);
  size_t taken = 0;
  for (size_t x = 0; x < redo.size(); ++x)
    if (!only_certified || res2[x].certified == 1) {
      S.res[redo[x]] = res2[x];
      ++taken;
    }
  return taken;
}

// The pair certificate on pair `pair` (reads 2 pair, 2 pair + 1): when it
// certifies, both results are replaced.
void MapPipeline::pair_certificate(Slot& S, size_t pair, Workspace& w1, Workspace& w2, PairScratch& sc) const {
  const ReadBatch& b = S.batch;
  const size_t i = 2 * pair, j = 2 * pair + 1;
  Result r1 = S.res[i], r2 = S.res[j];
  const PairOut po = certify_pair(view_, c_.pp, b.codes.data() + b.offs[i], b.lens[i], b.codes.data() + b.offs[j],
                                  b.lens[j], w1, w2, sc, r1, r2);
  if (po.certified) {
    S.res[i] = r1;
    S.res[j] = r2;
  }
}

RepeatOut MapPipeline::tier_r(const Slot& S, size_t i, Workspace& ws, RepeatScratch& sc) const {
  return repeat_tier(view_, c_.rp, S.batch.codes.data() + S.batch.offs[i], S.batch.lens[i], ws, sc);
}

// Formats reads [a0, b0) of a batch (whole pairs in paired mode: a0 and b0
// even) into `out`, in input order, and counts them into s; psc is scratch
// for interval mode's mate rescue.
void MapPipeline::format(const Slot& S, size_t a0, size_t b0, Chunk& out, Stats& s, PairScratch& psc) const {
  const std::vector<FastqRecord>& recs = S.recs;
  const std::vector<Result>& res = S.res;
  if (!c_.paired) {
    for (size_t i = a0; i < b0; ++i) {
      int reason;
      if (c_.intervals) {
        Placement pl;
        writer_.single_interval(res[i], S.iv_ptr[i], recs[i], out.sam, &pl);
        Result r = res[i];
        check_contig(r, ref_);
        reason = interval_reason(r);
        s.count_iv(pl);
      } else {
        reason = writer_.single(res[i], recs[i], out.sam, out.fq);
      }
      s.count_read(res[i], recs[i].seq.size(), reason);
    }
    return;
  }
  for (size_t i = a0; i < b0; i += 2) {
    int reason[2];
    ++s.pairs;
    if (c_.intervals) {
      const Result rr[2] = {mate_as_reported(res[i], ref_), mate_as_reported(res[i + 1], ref_)};
      Placement pl[2];
      if (rr[0].certified == 1 && rr[1].certified == 1) {
        writer_.pair(&res[i], &S.se_res[i], &recs[i], out.sam, out.fq, out.fq2, nullptr, reason, pl);
        ++s.pairs_cert;
        const bool proper = is_proper(rr[0], rr[1], ref_, c_.pp.max_dist);
        s.pairs_proper += proper;
        s.iv_proper += proper;
      } else {
        const RepeatOut* ivp[2] = {S.iv_ptr[i], S.iv_ptr[i + 1]};
        const uint8_t* cp[2] = {S.batch.codes.data() + S.batch.offs[i], S.batch.codes.data() + S.batch.offs[i + 1]};
        s.iv_proper += writer_.pair_interval(&res[i], ivp, cp, &recs[i], psc, out.sam, pl);
        for (int x = 0; x < 2; ++x) reason[x] = interval_reason(rr[x]);
      }
      s.count_iv(pl[0]);
      s.count_iv(pl[1]);
    } else {
      writer_.pair(&res[i], &S.se_res[i], &recs[i], out.sam, out.fq, out.fq2, out_.single ? &out.single : nullptr,
                   reason);
      if (reason[0] == kOk && reason[1] == kOk) {
        ++s.pairs_cert;
        s.pairs_proper += is_proper(res[i], res[i + 1], ref_, c_.pp.max_dist);
      }
    }
    for (int x = 0; x < 2; ++x) s.count_read(res[i + x], recs[i + x].seq.size(), reason[x]);
  }
}

// Hands the chunks to the writers, in order.
void MapPipeline::write(std::vector<Chunk>& chunks) {
  for (Chunk& ch : chunks) {
    out_.sam->put(std::move(ch.sam));
    if (out_.fq) out_.fq->put(std::move(ch.fq));
    if (out_.fq2) out_.fq2->put(std::move(ch.fq2));
    if (out_.single) out_.single->put(std::move(ch.single));
  }
}

// Passes 1 and 2 (on the GPU when there is one).
void MapPipeline::gpu_stage(Slot& S) {
  auto t = Clock::now();
  std::vector<Result>& res = S.res;
  if (gm_) gm_->map(c_.p1, S.batch, res);
  else map_cpu(view_, c_.p1, S.batch, res, c_.threads);
  if (c_.pass2()) {
    // Pass 2: every read pass 1 could not certify. Certified reads keep
    // their pass-1 results.
    std::vector<size_t> redo;
    for (size_t i = 0; i < res.size(); ++i)
      if (retry(res[i])) redo.push_back(i);
    rerun(S, redo, c_.p2, gm_ != nullptr);
    gpu_counts_.pass2 += redo.size();
  }
  time_.map_gpu += secs(t);
}

// The CPU passes, each over the whole batch, then formatting and writing, in
// input order (--serial-stages).
void MapPipeline::cpu_stage(Slot& S) {
  auto t = Clock::now();
  std::vector<Result>& res = S.res;
  PassCounts& n = cpu_counts_;
  auto reads_where = [&](auto want) {
    std::vector<size_t> idx;
    for (size_t i = 0; i < res.size(); ++i)
      if (want(res[i])) idx.push_back(i);
    return idx;
  };
  if (c_.pass3()) {
    // Pass 3, always on the CPU (identical in CPU and GPU runs).
    const std::vector<size_t> redo = reads_where([&](const Result& r) { return retry_pass3(r, c_.p3); });
    rerun(S, redo, c_.p3, false);
    n.pass3 += redo.size();
  }
  if (c_.pass_v()) {
    // Pass V: only a certificate replaces the pass-3 result.
    const std::vector<size_t> redo = reads_where(retry);
    n.pass_v_certified += rerun(S, redo, c_.pv, false, true);
    n.pass_v += redo.size();
  }
  if (c_.pass_l()) {
    // Pass L: certified local alignment (tier SL).
    const std::vector<size_t> redo = reads_where(retry_local);
    rerun(S, redo, c_.pl, false);
    n.pass_l += redo.size();
  }
  if (c_.paired) {
    S.se_res = res;
    std::vector<size_t> todo;  // pairs
    for (size_t i = 0; 2 * i + 1 < res.size(); ++i)
      if (wants_pair_pass(res[2 * i], res[2 * i + 1])) todo.push_back(i);
    n.pairs += todo.size();
    std::atomic<size_t> next{0};
    run_workers(c_.threads, [&](int) {
      std::unique_ptr<Workspace> w1(new Workspace), w2(new Workspace);
      PairScratch sc;
      for (size_t z; (z = next.fetch_add(1)) < todo.size();) pair_certificate(S, todo[z], *w1, *w2, sc);
    });
  }
  if (c_.intervals) {
    // Tier R for every read not certified.
    const std::vector<size_t> todo = reads_where([&](const Result& r) { return wants_tier_r(r, ref_); });
    S.iv_store.assign(todo.size(), RepeatOut());
    n.tier_r += todo.size();
    std::atomic<size_t> next{0};
    run_workers(c_.threads, [&](int) {
      std::unique_ptr<Workspace> w(new Workspace);
      RepeatScratch sc;
      for (size_t z; (z = next.fetch_add(1)) < todo.size();) S.iv_store[z] = tier_r(S, todo[z], *w, sc);
    });
    S.iv_ptr.assign(res.size(), nullptr);
    for (size_t z = 0; z < todo.size(); ++z) S.iv_ptr[todo[z]] = &S.iv_store[z];
  }
  time_.map_cpu += secs(t);

  // Format in parallel over contiguous chunks, which the threads take in
  // turn (hard reads cluster, so one fixed slice per thread leaves threads
  // idle), then write the chunks in input order.
  t = Clock::now();
  const size_t nr = S.recs.size();
  const int nt = std::max(1, std::min<int>(c_.threads, static_cast<int>(nr / 1024) + 1));
  const int nc = nt == 1 ? 1 : 8 * nt;
  std::vector<Chunk> chunks(nc);
  std::vector<Stats> st(nt);
  std::atomic<int> next_chunk{0};
  run_workers(nt, [&](int w) {
    PairScratch psc;  // interval mode: mate rescue
    for (int k; (k = next_chunk.fetch_add(1)) < nc;) {
      size_t a0 = nr * k / nc, b0 = nr * (k + 1) / nc;
      if (c_.paired) {
        a0 &= ~size_t(1);
        b0 = k == nc - 1 ? nr : (b0 & ~size_t(1));
      }
      format(S, a0, b0, chunks[k], st[w], psc);
    }
  });
  write(chunks);
  for (const Stats& s : st) total_.add(s);
  time_.output += secs(t);
}

// Every stage of a batch in one pass over chunks of reads (whole pairs),
// which the threads take in turn. Each chunk gets passes 1 and 2 (with
// `map`; with a GPU they ran in gpu_stage), passes 3, V and L, the second
// index, the pair pass, tier R and formatting, each read exactly what the
// staged passes give it (no result depends on which reads share a pass, a
// thread or a workspace). The threads then wait for each other once per
// batch instead of once per stage, so the slow reads of a stage no longer
// leave the other cores idle at its end.
void MapPipeline::fused_stage(Slot& S, bool map) {
  auto t = Clock::now();
  const size_t n = S.recs.size();
  S.res.resize(n);
  if (c_.paired) S.se_res.resize(n);
  if (c_.intervals) S.iv_ptr.assign(n, nullptr);
  const size_t nc = (n + kChunk - 1) / kChunk;
  std::vector<Chunk> chunks(nc);
  const int nt = std::max(1, std::min<int>(c_.threads, static_cast<int>(nc)));
  std::vector<Stats> st(nt);
  std::vector<PassCounts> counts(nt);
  std::vector<std::exception_ptr> errors(nt);
  std::atomic<size_t> next{0};
  run_workers(nt, [&](int w) {
    try {
      ChunkScratch x(view2_ != nullptr);
      for (size_t k; (k = next.fetch_add(1)) < nc;)
        fused_chunk(S, k * kChunk, std::min(n, (k + 1) * kChunk), map, x, chunks[k], st[w], counts[w]);
    } catch (...) {
      errors[w] = std::current_exception();
      next = nc;  // the others stop after their chunk
    }
  });
  for (auto& e : errors)
    if (e) std::rethrow_exception(e);
  time_.map_cpu += secs(t);
  t = Clock::now();
  write(chunks);
  for (int w = 0; w < nt; ++w) {
    total_.add(st[w]);
    cpu_counts_.add(counts[w]);
  }
  time_.output += secs(t);
}

// The fused stage on reads [a0, b0) of a batch.
void MapPipeline::fused_chunk(Slot& S, size_t a0, size_t b0, bool map, ChunkScratch& x, Chunk& out, Stats& s,
                              PassCounts& n) const {
  const ReadBatch& batch = S.batch;
  std::vector<Result>& res = S.res;
  size_t* const idx = x.idx;
  Result* const rr = x.rr.get();
  // idx[0..m) = the reads of the chunk that `want` selects; returns m.
  auto select = [&](auto want) {
    size_t m = 0;
    for (size_t i = a0; i < b0; ++i)
      if (want(res[i])) idx[m++] = i;
    return m;
  };
  // Re-runs reads idx[0..m) with p.
  auto rerun_chunk = [&](size_t m, const Params& p) {
    map_reads(view_, p, batch, idx, m, rr, *x.ws);
    for (size_t k = 0; k < m; ++k) res[idx[k]] = rr[k];
  };
  if (map) {
    size_t m = 0;
    for (size_t i = a0; i < b0; ++i) idx[m++] = i;
    map_reads(view_, c_.p1, batch, idx, m, res.data() + a0, *x.ws);
    if (c_.pass2()) {
      m = select(retry);
      rerun_chunk(m, c_.p2);
      n.pass2 += m;
    }
  }
  if (c_.pass3()) {
    const size_t m = select([&](const Result& r) { return retry_pass3(r, c_.p3); });
    rerun_chunk(m, c_.p3);
    n.pass3 += m;
  }
  if (c_.pass_v()) {
    // Pass V: a certificate found with adaptive parts is taken; the other
    // reads keep their pass-3 results (pass L reads those).
    const size_t m = select(retry);
    map_reads(view_, c_.pv, batch, idx, m, rr, *x.ws);
    for (size_t k = 0; k < m; ++k)
      if (rr[k].certified == 1) {
        res[idx[k]] = rr[k];
        ++n.pass_v_certified;
      }
    n.pass_v += m;
  }
  if (c_.pass_l()) {
    const size_t m = select(retry_local);
    rerun_chunk(m, c_.pl);
    n.pass_l += m;
  }
  if (view2_) index2_passes(S, a0, b0, x, n);
  if (c_.paired) {
    for (size_t i = a0; i < b0; ++i) S.se_res[i] = res[i];
    for (size_t i = a0 / 2; 2 * i + 1 < b0; ++i) {
      if (!wants_pair_pass(res[2 * i], res[2 * i + 1])) continue;
      ++n.pairs;
      pair_certificate(S, i, *x.w1, *x.w2, x.sc);
    }
  }
  if (c_.intervals) {
    x.civ.clear();
    const size_t m = select([&](const Result& r) { return wants_tier_r(r, ref_); });
    for (size_t k = 0; k < m; ++k) x.civ.push_back(tier_r(S, idx[k], *x.wr, x.rsc));
    for (size_t k = 0; k < m; ++k) S.iv_ptr[idx[k]] = &x.civ[k];
    n.tier_r += m;
  }
  format(S, a0, b0, out, s, x.psc);
}

// Every pass again on the second index (passes 1, 2, 3 and L) for the reads
// of [a0, b0) still uncertified; a certificate found there is taken.
void MapPipeline::index2_passes(Slot& S, size_t a0, size_t b0, ChunkScratch& x, PassCounts& n) const {
  const ReadBatch& batch = S.batch;
  std::vector<Result>& res = S.res;
  Result* const rx = x.rx.get();
  size_t m2 = 0;
  for (size_t i = a0; i < b0; ++i)
    if (retry(res[i])) x.id2[m2++] = i;
  if (m2 == 0) return;
  map_reads(*view2_, c_.p1, batch, x.id2, m2, rx, *x.wx);
  // One pass re-runs those of rx[0..m2) that `want` selects.
  auto rerun_index2 = [&](const Params& p, auto want) {
    size_t m = 0;
    for (size_t k = 0; k < m2; ++k)
      if (want(rx[k])) {
        x.sub[m] = k;
        x.idx[m++] = x.id2[k];
      }
    map_reads(*view2_, p, batch, x.idx, m, x.rr.get(), *x.wx);
    for (size_t j = 0; j < m; ++j) rx[x.sub[j]] = x.rr[j];
  };
  if (c_.pass2()) rerun_index2(c_.p2, retry);
  if (c_.pass3()) rerun_index2(c_.p3, [&](const Result& r) { return retry_pass3(r, c_.p3); });
  if (c_.pass_l()) rerun_index2(c_.pl, retry_local);
  for (size_t k = 0; k < m2; ++k)
    if (rx[k].certified == 1) {
      res[x.id2[k]] = rx[k];
      ++n.index2_certified;
    }
  n.index2 += m2;
}

void MapPipeline::run(const std::function<GpuMapper*()>& wait_for_gpu) {
  load(0);
  gm_ = wait_for_gpu();
  // Without a GPU each batch runs fused; with one, the CPU passes and
  // formatting do (overlapping the next batch's GPU passes).
  // --serial-stages runs the passes over the whole batch, one after the
  // other, instead.
  const bool fused = gm_ == nullptr && !c_.serial_stages;
  std::thread cpu_thread;
  std::exception_ptr cpu_error;
  Joiner cpu_joiner{cpu_thread};
  auto finish_cpu = [&] {
    if (cpu_thread.joinable()) cpu_thread.join();
    if (cpu_error) std::rethrow_exception(cpu_error);
  };
  for (int cur = 0; !slots_[cur].recs.empty(); cur = (cur + 1) % kSlots) {
    // The slot read into is neither the one the GPU stage maps nor the
    // one the CPU stage may still be writing (three slots).
    std::thread prefetch([this, cur] { load((cur + 1) % kSlots); });
    Joiner joiner{prefetch};
    Slot& S = slots_[cur];
    if (fused) {
      fused_stage(S, true);
    } else {
      gpu_stage(S);
      finish_cpu();  // the previous batch's CPU stage: batches are written in order
      if (gm_) {
        cpu_thread = std::thread([this, &S, &cpu_error] {
          try {
            if (c_.serial_stages) cpu_stage(S);
            else fused_stage(S, false);
          } catch (...) {
            cpu_error = std::current_exception();
          }
        });
      } else {
        cpu_stage(S);
      }
    }
    auto t = Clock::now();
    prefetch.join();
    time_.read_wait += secs(t);
  }
  finish_cpu();
  if (read_error_) std::rethrow_exception(read_error_);
}

}  // namespace certa
