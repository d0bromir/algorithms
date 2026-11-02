// The batches of certa map: read, map through every pass, write the records.
#pragma once
#include <cstddef>
#include <exception>
#include <functional>
#include <memory>
#include <vector>

#include "async_writer.h"
#include "certa/pair.h"
#include "certa/repeat.h"
#include "map_config.h"
#include "map_stats.h"
#include "mapper.h"
#include "sam_output.h"
#include "seqio.h"

namespace certa {

// The outputs (null: not written), each written by its own thread.
struct MapOutputs {
  std::unique_ptr<AsyncWriter> sam;
  std::unique_ptr<AsyncWriter> fq;      // reads for the fallback (mate 1)
  std::unique_ptr<AsyncWriter> fq2;     // their mates
  std::unique_ptr<AsyncWriter> single;  // --fallback-mates single: lone mates
  // Hands records to the writers, in order.
  void put(Records& r);
  // Writes and closes every output; false if a write failed.
  bool finish();
};

// Batches in flight: a reader thread parses and encodes the next batch while
// the current one maps. The first stage (passes 1 and 2) maps a batch over
// its whole length, on the GPU when there is one; the second (the CPU
// passes, formatting and writing) then finishes it, on its own thread with a
// GPU, so the GPU maps batch n+1 meanwhile. One second stage runs at a time,
// in input order, so the output does not depend on the overlap. Without a
// GPU (and without --serial-stages) every stage of a batch runs in one pass
// over chunks of reads instead (fused_stage).
class MapPipeline {
 public:
  // `view2`: the second index (--index2), or null; `mates`: null for
  // single-end input.
  MapPipeline(const MapConfig& c, const Reference& ref, const IndexView& view, const IndexView* view2,
              FastqInput& reads, FastqInput* mates, MapOutputs& out);

  // Maps and writes every batch. wait_for_gpu() is called once the first
  // batch is read (so reading overlaps the GPU upload) and returns the GPU
  // mapper, or null to map on the CPU.
  void run(const std::function<GpuMapper*()>& wait_for_gpu);

  const Stats& stats() const { return total_; }
  PassCounts counts() const;
  // Reading, mapping and output; load, upload and total are left to the caller.
  const Timings& timings() const { return time_; }

 private:
  // One batch: its records, their codes and results.
  struct Slot {
    std::vector<FastqRecord> recs;
    ReadBatch batch;
    std::vector<Result> res;
    std::vector<Result> se_res;            // single-end results (paired mode, for MAPQ)
    std::vector<RepeatOut> iv_store;       // tier R outputs (interval mode, serial stages)
    std::vector<const RepeatOut*> iv_ptr;  // per read: its tier R output, or null
    std::vector<BioOut> bio;               // --bio: per read, its claims
  };
  struct ChunkScratch;
  static constexpr int kSlots = 3;
  static constexpr size_t kChunk = 64;  // reads per fused chunk; even, so a pair never straddles two

  void load(int slot);
  size_t rerun(Slot& S, const std::vector<size_t>& redo, const Params& p, bool on_gpu, bool only_certified = false);
  void pair_certificate(Slot& S, size_t pair, Workspace& w1, Workspace& w2, PairScratch& sc) const;
  RepeatOut tier_r(const Slot& S, size_t i, Workspace& ws, RepeatScratch& sc) const;
  void format(const Slot& S, size_t a0, size_t b0, Records& out, Stats& s, PairScratch& psc) const;
  void write(std::vector<Records>& chunks);
  void first_stage(Slot& S);
  void cpu_stage(Slot& S);
  void fused_stage(Slot& S, bool map);
  void fused_chunk(Slot& S, size_t a0, size_t b0, bool map, ChunkScratch& x, Records& out, Stats& s,
                   PassCounts& n) const;
  void index2_passes(Slot& S, size_t a0, size_t b0, ChunkScratch& x, PassCounts& n) const;
  void bio_chunk(Slot& S, size_t a0, size_t b0, ChunkScratch& x) const;
  void bio_rescue_pair(Slot& S, size_t i, BioScratch& sc) const;

  const MapConfig& c_;
  const Reference& ref_;
  const IndexView& view_;
  const IndexView* view2_;
  FastqInput& reads_;
  FastqInput* mates_;
  MapOutputs& out_;
  RecordWriter writer_;
  BioModel bio_model_;  // --bio
  std::vector<uint8_t> bio_runs_;  // --bio: the reference's tandem runs (bio_pack_contexts), or empty
  GpuMapper* gm_ = nullptr;
  Slot slots_[kSlots];
  FastqSpans spans_[2];  // the mate files' text of the batch being read (buffers reused)
  std::exception_ptr read_error_;
  Stats total_;
  // Kept apart because with a GPU the two stages run on different threads.
  PassCounts first_counts_, rest_counts_;
  Timings time_;
};

}  // namespace certa
