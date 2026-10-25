// certa map: certified reads are written to SAM, all other reads unchanged to
// the uncertified FASTQ, for a full aligner (minibwa, BWA-MEM2, ...). With
// --intervals every read is written to SAM (certa/interval.h): the best
// alignment found, the interval [XL, XU] that contains its optimum, and XQ,
// a certified lower bound of its MAPQ.
#include <cstdio>
#include <cstring>
#include <future>
#include <memory>
#include <stdexcept>
#include <string>

#include "async_writer.h"
#include "cli.h"
#include "index.h"
#include "map_config.h"
#include "map_pipeline.h"
#include "map_stats.h"
#include "mapper.h"

namespace certa {

int cmd_map(const Args& a) {
  MapConfig c = MapConfig::from_args(a);

  auto t_all = Clock::now();
  Reference ref;
  Index ix;
  load_index(c.index_path, ref, ix);
  const double t_load = secs(t_all);
  std::fprintf(stderr, "[map] index %s: q=%d s=%d, %llu entries, loaded in %.1f s\n", c.index_path.c_str(), ix.q,
               ix.s, static_cast<unsigned long long>(ix.pos.size()), t_load);
  std::fprintf(stderr, "[map] pass 1: k=%d budget=%d; pass 2: k=%d budget=%d S2<=%d; parts are >= %d bases\n",
               c.p1.k, c.p1.budget, c.p2.k, c.p2.budget, c.p2.s2_limit, min_read_length(ix.q, ix.s));
  if (c.pass3() || c.pass_l())
    std::fprintf(stderr, "[map] pass 3: budget=%d support=%d; pass L: support %d%s\n", c.p3.budget,
                 c.p3.min_support, c.pl.local, c.full_cert ? " (--full-cert)" : "");
  c.fit_index(ix.s);
  if (c.pass_v())
    std::fprintf(stderr, "[map] pass V: adaptive parts %s, budget=%d, radius <= %d\n",
                 c.pv.vparts == 1 ? "(rarest q-mers)" : "(rarest q-mers, then one-edit q-mers)", c.pv.budget, VMAX);
  const IndexView view = ix.view(ref);
  // The second index must hold the same reference: its view reads the first
  // index's copy of the sequence.
  Index ix2;
  std::unique_ptr<IndexView> view2;
  if (c.index2()) {
    auto t = Clock::now();
    Reference ref2;
    load_index(c.index2_path, ref2, ix2);
    if (ref2.names != ref.names || ref2.lengths != ref.lengths || ref2.seq.size() != ref.seq.size() ||
        std::memcmp(ref2.seq.data(), ref.seq.data(), ref.seq.size()) != 0)
      throw std::runtime_error("--index2 " + c.index2_path + " was built from another reference");
    view2.reset(new IndexView(ix2.view(ref)));
    std::fprintf(stderr, "[map] second index %s: q=%d s=%d, %llu entries, loaded and checked in %.1f s\n",
                 c.index2_path.c_str(), ix2.q, ix2.s, static_cast<unsigned long long>(ix2.pos.size()), secs(t));
  }

  // Open input and outputs first, so a failure cannot leave the upload
  // thread running.
  FastqInput reads(c.reads_path);
  std::unique_ptr<FastqInput> mates(c.paired() ? new FastqInput(c.mates_path) : nullptr);
  std::FILE* sam = c.sam_path.empty() ? stdout : open_output(c.sam_path, "SAM output");
  std::FILE* fq = c.fq_path.empty() ? nullptr : open_output(c.fq_path, "uncertified FASTQ output");
  std::FILE* fq2 = c.fq2_path.empty() ? nullptr : open_output(c.fq2_path, "uncertified mates output");
  std::FILE* fqs = c.single_path.empty() ? nullptr : open_output(c.single_path, "single-mate fallback output");
  std::fprintf(sam, "@HD\tVN:1.6\tSO:unsorted\n");
  for (size_t i = 0; i < ref.names.size(); ++i)
    std::fprintf(sam, "@SQ\tSN:%s\tLN:%llu\n", ref.names[i].c_str(), static_cast<unsigned long long>(ref.lengths[i]));
  std::fprintf(sam, "@PG\tID:certa\tPN:certa\tVN:0.5\tCL:%s\n", c.command_line().c_str());
  // From here on each output is written (and closed) by its own thread. Two
  // mate files read in lockstep by one consumer cannot deadlock: each writer
  // moves on as soon as its own pipe has room.
  constexpr size_t kWriteQueue = size_t(2) << 30;
  auto writer = [](std::FILE* f) { return std::unique_ptr<AsyncWriter>(f ? new AsyncWriter(f, kWriteQueue) : nullptr); };
  MapOutputs out{writer(sam), writer(fq), writer(fq2), writer(fqs)};

  // The GPU upload runs in the background, overlapped with reading the
  // first batch; the first map call waits for it.
  std::unique_ptr<GpuMapper> gm;
  double t_upload = 0;
  std::future<void> upload;
  if (c.gpu)
    upload = std::async(std::launch::async, [&] {
      auto t = Clock::now();
      gm.reset(new GpuMapper(ref, ix, c.device));
      t_upload = secs(t);
    });
  auto wait_for_gpu = [&]() -> GpuMapper* {
    if (!upload.valid()) return nullptr;
    upload.get();  // rethrows an upload error
    std::fprintf(stderr, "[map] GPU %s, index uploaded in %.1f s (overlapped with input)\n",
                 gm->device_name().c_str(), t_upload);
    return gm.get();
  };

  MapPipeline pipeline(c, ref, view, view2.get(), reads, mates.get(), out);
  pipeline.run(wait_for_gpu);
  if (!out.finish()) throw std::runtime_error("writing an output failed");

  Timings t = pipeline.timings();
  t.load = t_load;
  t.upload = t_upload;
  t.total = secs(t_all);
  print_summary(c, pipeline.stats(), pipeline.counts(), t);
  if (!c.stats_path.empty())
    write_stats_json(c, ix, gm ? gm->device_name() : "", pipeline.stats(), pipeline.counts(), t);
  return 0;
}

}  // namespace certa
