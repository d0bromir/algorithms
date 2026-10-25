// The batches of certa map: read, map through every pass, write the records.
#pragma once
#include <cstddef>
#include <exception>
#include <functional>
#include <memory>
#include <string>
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

// A FASTQ input: memory-mapped when it is an uncompressed regular file,
// streamed otherwise (gzip, pipes).
class FastqInput {
 public:
  explicit FastqInput(const std::string& path) {
    if (MappedFastq::usable(path)) mapped_.reset(new MappedFastq(path));
    else stream_.reset(new LineReader(path));
  }
  // Up to n records into `out`; they are built on `threads` threads.
  size_t next_batch(std::vector<FastqRecord>& out, size_t n, int threads) {
    return mapped_ ? mapped_->next_batch(out, n, threads) : read_fastq_batch(*stream_, out, n, threads);
  }
  // Up to n records as text, for parse_spans (a stream is split on one thread).
  size_t next_spans(FastqSpans& sp, size_t n, int threads) {
    return mapped_ ? mapped_->next_spans(sp, n, threads) : read_fastq_spans(*stream_, sp, n);
  }

 private:
  std::unique_ptr<MappedFastq> mapped_;
  std::unique_ptr<LineReader> stream_;
};

// The outputs (null: not written).
struct MapOutputs {
  AsyncWriter* sam = nullptr;
  AsyncWriter* fq = nullptr;      // reads for the fallback (mate 1)
  AsyncWriter* fq2 = nullptr;     // their mates
  AsyncWriter* single = nullptr;  // --fallback-mates single: lone mates
};

// Batches in flight: a reader thread parses and encodes the next batch while
// the current one maps. With a GPU, the first stage (passes 1 and 2, on the
// GPU) maps batch n+1 while the second (the CPU passes, formatting and
// writing) finishes batch n; one second stage runs at a time, in input
// order, so the output does not depend on the overlap. Without a GPU every
// stage of a batch runs in one pass over chunks of reads (fused_stage).
class MapPipeline {
 public:
  // `view2`: the second index (--index2), or null; `mates`: null for
  // single-end input.
  MapPipeline(const MapConfig& c, const Reference& ref, const IndexView& view, const IndexView* view2,
              FastqInput& reads, FastqInput* mates, const MapOutputs& out);

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
  };
  // The records of a range of reads, one string per output.
  struct Chunk {
    std::string sam, fq, fq2, single;
  };
  struct ChunkScratch;
  static constexpr int kSlots = 3;
  static constexpr size_t kChunk = 64;  // reads per fused chunk; even, so a pair never straddles two

  void load(int slot);
  size_t rerun(Slot& S, const std::vector<size_t>& redo, const Params& p, bool on_gpu, bool only_certified = false);
  void pair_certificate(Slot& S, size_t pair, Workspace& w1, Workspace& w2, PairScratch& sc) const;
  RepeatOut tier_r(const Slot& S, size_t i, Workspace& ws, RepeatScratch& sc) const;
  void format(const Slot& S, size_t a0, size_t b0, Chunk& out, Stats& s, PairScratch& psc) const;
  void write(std::vector<Chunk>& chunks);
  void gpu_stage(Slot& S);
  void cpu_stage(Slot& S);
  void fused_stage(Slot& S, bool map);
  void fused_chunk(Slot& S, size_t a0, size_t b0, bool map, ChunkScratch& x, Chunk& out, Stats& s,
                   PassCounts& n) const;
  void index2_passes(Slot& S, size_t a0, size_t b0, ChunkScratch& x, PassCounts& n) const;

  const MapConfig& c_;
  const Reference& ref_;
  const IndexView& view_;
  const IndexView* view2_;
  FastqInput& reads_;
  FastqInput* mates_;
  MapOutputs out_;
  RecordWriter writer_;
  GpuMapper* gm_ = nullptr;
  Slot slots_[kSlots];
  std::exception_ptr read_error_;
  Stats total_;
  // Kept apart because with a GPU the two stages run on different threads.
  PassCounts gpu_counts_, cpu_counts_;
  Timings time_;
};

}  // namespace certa
