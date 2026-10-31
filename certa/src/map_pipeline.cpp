#include "map_pipeline.h"

#include <algorithm>
#include <atomic>
#include <future>
#include <numeric>
#include <stdexcept>

#include "cli.h"

namespace certa {

namespace {

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

// The reads of a whole batch that `want` selects.
template <class Want>
std::vector<size_t> reads_where(const std::vector<Result>& res, Want want) {
  std::vector<size_t> idx;
  for (size_t i = 0; i < res.size(); ++i)
    if (want(res[i])) idx.push_back(i);
  return idx;
}

// The reads of [a0, b0) that `want` selects, into idx; returns how many.
template <class Want>
size_t select(const std::vector<Result>& res, size_t a0, size_t b0, Want want, size_t* idx) {
  size_t m = 0;
  for (size_t i = a0; i < b0; ++i)
    if (want(res[i])) idx[m++] = i;
  return m;
}

}  // namespace

void MapOutputs::put(Records& r) {
  sam->put(std::move(r.sam));
  if (fq) fq->put(std::move(r.fq));
  if (fq2) fq2->put(std::move(r.fq2));
  if (single) single->put(std::move(r.single));
}

bool MapOutputs::finish() {
  bool ok = true;
  for (auto* w : {sam.get(), fq.get(), fq2.get(), single.get()})
    if (w) ok = w->finish() && ok;
  return ok;
}

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
  std::unique_ptr<BioScratch> bsc{new BioScratch};  // --bio
};

MapPipeline::MapPipeline(const MapConfig& c, const Reference& ref, const IndexView& view, const IndexView* view2,
                         FastqInput& reads, FastqInput* mates, MapOutputs& out)
    : c_(c), ref_(ref), view_(view), view2_(view2), reads_(reads), mates_(mates), out_(out),
      writer_(ref, view, c.out), bio_model_(c.bp) {}

PassCounts MapPipeline::counts() const {
  PassCounts n = first_counts_;
  n.add(rest_counts_);
  return n;
}

// Reads and encodes the next batch into a slot (empty at the end of input).
// An error is kept for run() to throw.
void MapPipeline::load(int slot) {
  Slot& S = slots_[slot];
  try {
    auto t = Clock::now();
    if (!c_.paired()) {
      reads_.next_batch(S.recs, c_.batch_size, c_.io_threads);
    } else {
      // Mates interleaved: read 2i is mate 1, read 2i + 1 its mate 2, each
      // parsed straight into its place (strings reused from the slot's
      // previous batch). The two files are read at the same time (each is
      // split, and a gzip file inflated, on one thread).
      const size_t half = std::max<size_t>(1, c_.batch_size / 2);
      auto n2 = std::async(std::launch::async, [&] { return mates_->next_spans(spans_[1], half, c_.io_threads); });
      const size_t n1 = reads_.next_spans(spans_[0], half, c_.io_threads);
      if (n1 != n2.get()) throw std::runtime_error("mate files have different numbers of reads");
      std::vector<FastqRecord>& out = S.recs;
      out.resize(2 * n1);
      parse_spans(spans_[0], out.data(), 2, c_.io_threads);
      parse_spans(spans_[1], out.data() + 1, 2, c_.io_threads);
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
// certifies, both results are replaced, and the single-end ones kept in
// se_res for MAPQ (only tier PR reads them).
void MapPipeline::pair_certificate(Slot& S, size_t pair, Workspace& w1, Workspace& w2, PairScratch& sc) const {
  const ReadBatch& b = S.batch;
  const size_t i = 2 * pair, j = 2 * pair + 1;
  Result r1 = S.res[i], r2 = S.res[j];
  const PairOut po = certify_pair(view_, c_.pp, b.codes.data() + b.offs[i], b.lens[i], b.codes.data() + b.offs[j],
                                  b.lens[j], w1, w2, sc, r1, r2);
  if (po.certified) {
    S.se_res[i] = S.res[i];
    S.se_res[j] = S.res[j];
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
void MapPipeline::format(const Slot& S, size_t a0, size_t b0, Records& out, Stats& s, PairScratch& psc) const {
  const std::vector<FastqRecord>& recs = S.recs;
  const std::vector<Result>& res = S.res;
  if (c_.bio) {
    for (size_t i = a0; i < b0; i += c_.paired() ? 2 : 1) {
      if (!c_.paired()) {
        const int reason = writer_.bio_single(res[i], S.bio[i], recs[i], out);
        s.count_read(res[i], recs[i].seq.size(), reason);
        s.count_bio(res[i], S.bio[i], reason);
        continue;
      }
      int reason[2];
      const bool proper = writer_.bio_pair(&res[i], &S.bio[i], &recs[i], out, reason);
      ++s.pairs;
      if (reason[0] == kOk && reason[1] == kOk) {
        ++s.pairs_cert;
        s.pairs_proper += proper;
      }
      for (int x = 0; x < 2; ++x) {
        s.count_read(res[i + x], recs[i + x].seq.size(), reason[x]);
        s.count_bio(res[i + x], S.bio[i + x], reason[x]);
      }
    }
    return;
  }
  if (!c_.paired()) {
    for (size_t i = a0; i < b0; ++i) {
      int reason;
      if (c_.intervals) {
        Placement pl;
        reason = writer_.single_interval(res[i], S.iv_ptr[i], recs[i], out, &pl);
        s.count_iv(pl);
      } else {
        reason = writer_.single(res[i], recs[i], out);
      }
      s.count_read(res[i], recs[i].seq.size(), reason);
    }
    return;
  }
  for (size_t i = a0; i < b0; i += 2) {
    int reason[2];
    bool proper;
    if (c_.intervals) {
      Placement pl[2];
      const RepeatOut* ivp[2] = {S.iv_ptr[i], S.iv_ptr[i + 1]};
      const uint8_t* cp[2] = {S.batch.codes.data() + S.batch.offs[i], S.batch.codes.data() + S.batch.offs[i + 1]};
      proper = writer_.pair_interval(&res[i], &S.se_res[i], ivp, cp, &recs[i], psc, out, reason, pl);
      s.iv_proper += proper;
      s.count_iv(pl[0]);
      s.count_iv(pl[1]);
    } else {
      proper = writer_.pair(&res[i], &S.se_res[i], &recs[i], out, reason);
    }
    ++s.pairs;
    if (reason[0] == kOk && reason[1] == kOk) {
      ++s.pairs_cert;
      s.pairs_proper += proper;
    }
    for (int x = 0; x < 2; ++x) s.count_read(res[i + x], recs[i + x].seq.size(), reason[x]);
  }
}

void MapPipeline::write(std::vector<Records>& chunks) {
  for (Records& r : chunks) out_.put(r);
}

// Passes 1 and 2 over the whole batch (on the GPU when there is one).
void MapPipeline::first_stage(Slot& S) {
  auto t = Clock::now();
  std::vector<Result>& res = S.res;
  if (gm_) gm_->map(c_.p1, S.batch, res);
  else map_cpu(view_, c_.p1, S.batch, res, c_.threads);
  if (c_.pass2()) {
    // Pass 2: every read pass 1 could not certify. Certified reads keep
    // their pass-1 results.
    const std::vector<size_t> redo = reads_where(res, retry);
    rerun(S, redo, c_.p2, gm_ != nullptr);
    first_counts_.pass2 += redo.size();
  }
  time_.map_first += secs(t);
}

// The CPU passes, each over the whole batch, then formatting and writing, in
// input order (with a GPU, or --serial-stages).
void MapPipeline::cpu_stage(Slot& S) {
  auto t = Clock::now();
  std::vector<Result>& res = S.res;
  PassCounts& n = rest_counts_;
  if (c_.pass3()) {
    // Pass 3, always on the CPU (identical in CPU and GPU runs).
    const std::vector<size_t> redo = reads_where(res, [&](const Result& r) { return retry_pass3(r, c_.p3); });
    rerun(S, redo, c_.p3, false);
    n.pass3 += redo.size();
  }
  if (c_.pass_v()) {
    // Pass V: only a certificate replaces the pass-3 result.
    const std::vector<size_t> redo = reads_where(res, retry);
    n.pass_v_certified += rerun(S, redo, c_.pv, false, true);
    n.pass_v += redo.size();
  }
  if (c_.pass_l()) {
    // Pass L: certified local alignment (tier SL).
    const std::vector<size_t> redo = reads_where(res, retry_local);
    rerun(S, redo, c_.pl, false);
    n.pass_l += redo.size();
  }
  if (c_.paired()) {
    S.se_res.resize(res.size());
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
    const std::vector<size_t> todo = reads_where(res, [&](const Result& r) { return wants_tier_r(r, ref_); });
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
  time_.map_rest += secs(t);

  // Format in parallel over contiguous chunks, which the threads take in
  // turn (hard reads cluster, so one fixed slice per thread leaves threads
  // idle), then write the chunks in input order.
  t = Clock::now();
  const size_t nr = S.recs.size();
  const int nt = std::max(1, std::min<int>(c_.threads, static_cast<int>(nr / 1024) + 1));
  const int nc = nt == 1 ? 1 : 8 * nt;
  std::vector<Records> chunks(nc);
  std::vector<Stats> st(nt);
  std::atomic<int> next_chunk{0};
  run_workers(nt, [&](int w) {
    PairScratch psc;  // interval mode: mate rescue
    for (int k; (k = next_chunk.fetch_add(1)) < nc;) {
      size_t a0 = nr * k / nc, b0 = nr * (k + 1) / nc;
      if (c_.paired()) {
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
// `map`; with a GPU they ran in first_stage), passes 3, V and L, the second
// index, the pair pass, tier R and formatting, each read exactly what the
// staged passes give it (no result depends on which reads share a pass, a
// thread or a workspace). The threads then wait for each other once per
// batch instead of once per stage, so the slow reads of a stage no longer
// leave the other cores idle at its end.
void MapPipeline::fused_stage(Slot& S, bool map) {
  auto t = Clock::now();
  const size_t n = S.recs.size();
  S.res.resize(n);
  if (c_.paired()) S.se_res.resize(n);
  if (c_.intervals) S.iv_ptr.assign(n, nullptr);
  if (c_.bio) S.bio.assign(n, BioOut());
  const size_t nc = (n + kChunk - 1) / kChunk;
  std::vector<Records> chunks(nc);
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
  time_.map_rest += secs(t);
  t = Clock::now();
  write(chunks);
  for (int w = 0; w < nt; ++w) {
    total_.add(st[w]);
    rest_counts_.add(counts[w]);
  }
  time_.output += secs(t);
}

// The fused stage on reads [a0, b0) of a batch.
void MapPipeline::fused_chunk(Slot& S, size_t a0, size_t b0, bool map, ChunkScratch& x, Records& out, Stats& s,
                              PassCounts& n) const {
  const ReadBatch& batch = S.batch;
  std::vector<Result>& res = S.res;
  size_t* const idx = x.idx;
  Result* const rr = x.rr.get();
  // Re-runs reads idx[0..m) with p (with only_certified, keeping only
  // certificates); returns the results taken.
  auto rerun_chunk = [&](size_t m, const Params& p, bool only_certified = false) {
    map_reads(view_, p, batch, idx, m, rr, *x.ws);
    size_t taken = 0;
    for (size_t k = 0; k < m; ++k)
      if (!only_certified || rr[k].certified == 1) {
        res[idx[k]] = rr[k];
        ++taken;
      }
    return taken;
  };
  if (c_.bio) {  // one pass: the maximum-likelihood placement
    bio_chunk(S, a0, b0, x);
    format(S, a0, b0, out, s, x.psc);
    return;
  }
  if (map) {
    std::iota(idx, idx + (b0 - a0), a0);
    map_reads(view_, c_.p1, batch, idx, b0 - a0, res.data() + a0, *x.ws);
    if (c_.pass2()) {
      const size_t m = select(res, a0, b0, retry, idx);
      rerun_chunk(m, c_.p2);
      n.pass2 += m;
    }
  }
  if (c_.pass3()) {
    const size_t m = select(res, a0, b0, [&](const Result& r) { return retry_pass3(r, c_.p3); }, idx);
    rerun_chunk(m, c_.p3);
    n.pass3 += m;
  }
  if (c_.pass_v()) {
    // Pass V: a certificate found with adaptive parts is taken; the other
    // reads keep their pass-3 results (pass L reads those).
    const size_t m = select(res, a0, b0, retry, idx);
    n.pass_v_certified += rerun_chunk(m, c_.pv, true);
    n.pass_v += m;
  }
  if (c_.pass_l()) {
    const size_t m = select(res, a0, b0, retry_local, idx);
    rerun_chunk(m, c_.pl);
    n.pass_l += m;
  }
  if (view2_) index2_passes(S, a0, b0, x, n);
  if (c_.paired()) {
    for (size_t i = a0 / 2; 2 * i + 1 < b0; ++i) {
      if (!wants_pair_pass(res[2 * i], res[2 * i + 1])) continue;
      ++n.pairs;
      pair_certificate(S, i, *x.w1, *x.w2, x.sc);
    }
  }
  if (c_.intervals) {
    x.civ.clear();
    const size_t m = select(res, a0, b0, [&](const Result& r) { return wants_tier_r(r, ref_); }, idx);
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
  const size_t m2 = select(res, a0, b0, retry, x.id2);
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

// --bio: every read of [a0, b0) through bio_place. In paired mode, a pair
// whose mates both have a certified ML placement is then placed jointly
// when the pair's optimum is certified (bio_pair): each mate is re-aligned
// at its locus in that pair, and takes the pair's posterior as MAPQ.
void MapPipeline::bio_chunk(Slot& S, size_t a0, size_t b0, ChunkScratch& x) const {
  const ReadBatch& b = S.batch;
  BioScratch& sc = *x.bsc;
  for (size_t i = a0; i < b0; ++i)
    bio_place(view_, bio_model_, c_.pb, b.codes.data() + b.offs[i], S.recs[i].qual.data(), b.lens[i], b.hashes[i],
              *x.ws, sc, S.res[i], S.bio[i]);
  if (!c_.paired()) return;
  auto contig_of = [&](int64_t pos) { return ref_.contig_of(pos); };
  for (size_t i = a0; i + 1 < b0; i += 2) {
    BioOut* o = &S.bio[i];
    const int L[2] = {b.lens[i], b.lens[i + 1]};
    const BioPair bp = bio_pair(c_.bp, o[0], L[0], o[1], L[1], c_.pp.max_dist, contig_of, b.hashes[i]);
    if (!bp.ok) continue;
    const int pick[2] = {bp.x, bp.y};
    Result r[2] = {S.res[i], S.res[i + 1]};
    bool aligned = true;
    for (int k = 0; k < 2 && aligned; ++k) {
      if (pick[k] == o[k].chosen) continue;
      const size_t j = i + k;
      bio_read(bio_model_, b.codes.data() + b.offs[j], S.recs[j].qual.data(), L[k], sc.br);
      aligned = bio_align_locus(view_, bio_model_, o[k].loci[pick[k]], sc, r[k]);
    }
    if (!aligned) continue;
    for (int k = 0; k < 2; ++k) {
      S.res[i + k] = r[k];
      o[k].chosen = pick[k];
      o[k].pair_ok = true;
      o[k].pair_pen = bp.pen;
      o[k].ins_pen = bp.ins;
      o[k].pair_mapq = bp.mapq[k];
    }
  }
}

void MapPipeline::run(const std::function<GpuMapper*()>& wait_for_gpu) {
  load(0);
  gm_ = wait_for_gpu();
  // Passes 1 and 2 run over the whole batch first with a GPU or
  // --serial-stages; the rest then runs serially (--serial-stages) or fused,
  // on its own thread when a GPU maps the next batch meanwhile. Without
  // either, every pass runs fused.
  const bool split = gm_ != nullptr || c_.serial_stages;
  auto rest = [this, split](Slot& S) {
    if (c_.serial_stages) cpu_stage(S);
    else fused_stage(S, !split);
  };
  std::future<void> cpu;  // the previous batch's second stage (with a GPU)
  for (int cur = 0; !slots_[cur].recs.empty(); cur = (cur + 1) % kSlots) {
    // The slot read into is neither the one the first stage maps nor the
    // one the second stage may still be writing (three slots).
    auto prefetch = std::async(std::launch::async, [this, cur] { load((cur + 1) % kSlots); });
    Slot& S = slots_[cur];
    if (split) first_stage(S);
    if (cpu.valid()) cpu.get();  // batches are written in order
    if (gm_) cpu = std::async(std::launch::async, rest, std::ref(S));
    else rest(S);
    auto t = Clock::now();
    prefetch.wait();
    time_.read_wait += secs(t);
  }
  if (cpu.valid()) cpu.get();
  if (read_error_) std::rethrow_exception(read_error_);
}

}  // namespace certa
