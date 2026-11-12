# One pass from FASTQ to certified alignments: fastp-gpu × CERTA (design)

**Status:** design (2026-10-29). M0, the fastp library, is built and tested on
CPU (2026-11-07): branch `claude/m0-libfastp-core` of `d0bromir/fastp`, the
work repository that `fastp-gpu` is published from (see its `docs/LIBRARY.md`).
Sections 3.3 and 5 describe the library as built. A first M1 is built
(2026-11-12): `certa-qc` in `certa/`, built when CMake is given `FASTP_DIR`.
Its output equals the two-step pipeline. It still uses two thread pools and
writes reads hard-trimmed (section 3.4's soft clips are not done yet). On 4
cores it is about as fast as fastp, a plain file, then certa map (section 8).
Repositories:
[d0bromir/fastp-gpu](https://github.com/d0bromir/fastp-gpu) (v1.3.3-d0bromir,
QC/trimming, CUDA, nvCOMP, GPU-Direct Storage) and this repository's
[`certa/`](certa/README.md) (certified mapping, CUDA). Working name of the new
program: **`certa-qc`**. The product built on it (assay detection,
fit / not-fit reports, contamination) is planned in
[FIT_FOR_PURPOSE_PIPELINE.md](FIT_FOR_PURPOSE_PIPELINE.md).

## 1. Goal

One process that reads the FASTQ files once, runs fastp's QC, trimming and
filtering and CERTA's certified mapping on the same reads, and writes the
alignments, the fallback reads and both reports:

- **no intermediate FASTQ**: nothing is compressed, written, read back or
  inflated between the two programs;
- **one upload per batch**: a batch's text goes to the GPU once (or is
  produced there, from compressed input), QC and mapping read the same
  device buffer;
- **one CUDA context and one index upload**: the reference index (GRCh38:
  9–21 GiB) is uploaded once, overlapped with reading the first batch, as
  `certa map --gpu` already does;
- **the same output as the two programs run one after the other**, byte for
  byte, on any thread count, batch size and back-end (CPU or GPU). Both
  projects already test themselves this way; the fused program is tested
  against the pair of them.

## 2. What the hand-off costs today

The two-step pipeline today:

```
R1/R2.fq.gz ──inflate──▶ fastp-gpu ──deflate──▶ trimmed.fq.gz ──inflate──▶ certa map ──▶ certified.sam
                         (QC, trim,                             (parse, encode,          + uncertified.fq
                          filter)                                upload, map)              ──▶ fallback aligner
```

The input is inflated twice and the trimmed reads deflated once, and both
programs parse the text and upload it to the GPU separately. Measured in the
two repositories:

- **fastp-gpu is bound by inflation and output compression**
  ([docs/PERFORMANCE_SUMMARY.md](https://github.com/d0bromir/fastp-gpu/blob/main/docs/PERFORMANCE_SUMMARY.md)):
  inflating one 41.3 GB mate takes 204.5 s on one thread, and the CPU build
  reaches 241 s from 16 threads on, close to that bound; parallel output
  compression is where most of its 1.5–2.9× speedup over upstream comes from
  (2.4× at 8 threads, 4.4× at 32). The GPU QC kernels themselves are not the
  bottleneck (the GPU build is faster in 2 of 12 cells).
- **certa map** inflates a gzip input on one thread per file, and before the
  parallel reader its end-to-end time was dominated by parsing and encoding
  ([certa/README.md](certa/README.md), "Parsing now dominates").

So in the fused program the trimmed reads are never deflated or inflated,
the original input is inflated once, and the text is parsed once. With
uncompressed or BGZF input (parallel or GPU inflation) the remaining bound is
the mapping itself.

## 3. Architecture

```mermaid
flowchart LR
    subgraph IN["Input (once)"]
        A1["R1/R2 .fq.gz / BGZF / plain"]
        A2["ISA-L inflate (CPU)<br/>or nvCOMP (GPU)<br/>or GDS NVMe→GPU"]
        A1 --> A2
    end
    A2 --> B["Batch: text + record descriptors<br/>(host pinned or device),<br/>mates in lockstep"]
    B --> C["fastp per-read work<br/>stats, adapter/overlap, polyG/X,<br/>quality window, filter<br/>(GPU kernels + CPU)"]
    C --> D["Bridge: trim windows → 2-bit codes<br/>+ qualities (encode kernel,<br/>no host round trip)"]
    D --> E["CERTA passes 1–2 (GPU)<br/>3, V, L, pairs, tier R, --bio (CPU)"]
    E --> F["Writers: BAM/SAM (BGZF pool or nvCOMP),<br/>fallback FASTQ (~5 % of reads),<br/>fastp JSON/HTML + certa stats JSON"]
    C -- "failed reads" --> F
```

### 3.1 Batches

One batch type flows through every stage: the decompressed FASTQ text, a
record descriptor per read (`GpuReadDescriptor` in fastp-gpu: sequence offset,
length, quality offset), and per-read results filled in by each stage. In
paired mode both files are read in lockstep, so mate *i* of both files is in
the same batch (fastp's `PairEndProcessor` and `certa map` already require
this). Batches are numbered; writers emit them in order (fastp-gpu's
min-heap ordered writer, `certa`'s in-order `AsyncWriter`).

### 3.2 Stage 1: input

fastp-gpu's three readers are kept as they are: ISA-L inflation on the CPU,
nvCOMP BGZF/gzip inflation on the GPU (`cuda_gzip.cu`), and the GDS path
(`gds_pipeline.cu`: NVMe → GPU memory → BGZF parse → nvCOMP → FASTQ parse →
`GpuReadDescriptor[]`, all on the device). With the last two the batch text
is produced on the GPU and never crosses PCIe in the host-to-device
direction.

### 3.3 Stage 2: QC (fastp)

fastp's per-pack work runs unchanged: UMI, front/tail trimming, adapter
trimming by sequence or PE overlap, polyG/polyX, quality-window trimming,
filtering, deduplication, pre- and post-filter statistics. M0 moved it into
`SeCore` (single-end) and `PeCore` (paired-end). Where fastp appended a read
to an output string, the core now calls an output object, in the same order. fastp's
object appends the same text; `certa-qc`'s receives the reads.

**Where the trimming happens.** All trimming runs on the CPU and edits the
read in place. The GPU kernel's `trim_start`, `trim_end` and `polyG_trim_pos`
(`ReadStatistics`) are computed but not used to trim. So the trim window comes
from the CPU steps, not from the kernel.

**The trim window.** For each written read the core reports a front offset:
the bases removed from the read start by UMI extraction, `--trim_front` and
`--cut_front`. Every other step trims the end only. The written read equals
the input read from that offset on, for its length, except in two cases:

- `--correction` can change bases in the overlap. The core then marks both
  mates as modified.
- Merged reads (`-m`) and overlapped regions are new reads, always modified.

A modified read has to reach the mapper as text, not as a window. Because
fastp trims in place, the input read is gone after trimming. The core can keep
a copy of it on request, for the soft clips (section 3.4).

**The evaluator** (read length, adapter detection, two-colour check) reads
the start of the input files before processing, as in fastp. A library
session does the same, so its options equal fastp's.

### 3.4 Stage 3: the bridge

A kernel (or, without a GPU, the existing `encode_batch` loop) writes the
2-bit codes of each kept read's window straight into the mapper's code
buffer, with its offset, length and name hash, so `GpuMapper` maps a
device-resident batch. The windows come from fastp's CPU trimming (section
3.3) and are uploaded with the batch. Modified reads are uploaded as text.
The one new entry point in CERTA:

```cpp
// certa/src/mapper.h
struct DeviceReadBatch { const uint8_t* text; const GpuReadDescriptor* desc;
                         const int32_t* trim_start; const int32_t* trim_end;
                         const uint8_t* keep; size_t n; };
void GpuMapper::map_device(const Params& p, const DeviceReadBatch& b, std::vector<Result>& out);
```

Only the uncertified reads' codes and qualities come back to the host for the
CPU passes (about 5 % of reads by default, and the qualities of every read
with `--bio`), with the 200-byte `Result` per read.

**Trimmed bases in the output.** By default the SAM/BAM record keeps the
original read, and the trimmed bases are soft clips (CIGAR `S`) around the
alignment of the window: downstream tools see the whole read, and the
certificate is about the window, as CERTA's tags already say for the aligned
span. `--hard-trim` writes the trimmed read instead, as `fastp | certa`
would.

### 3.5 Stage 4: mapping

`certa map`'s pipeline (`MapPipeline`) is unchanged apart from where its
batches come from: an interface `BatchSource` replaces `FastqInput`. Passes 1
and 2 run on the GPU (or CPU), passes 3, V, L, the pair certificate, tier R
and `--bio` on the CPU, fused per chunk of reads as now.

Two things the fused program can do that two programs cannot:

- **The fragment-length model of `--bio` pairs from the data.** fastp's
  overlap analysis measures the insert size of every overlapping pair before
  any mapping. The first N pairs' histogram (N fixed, for reproducibility)
  replaces the `--bio-insert MEAN,SD` Gaussian with the empirical
  distribution.
- **Adapters and poly-G trimmed before certification.** `--bio` aligns end to
  end, so a read with an adapter or a poly-G tail costs too much and goes to
  the fallback. fastp trims both first; on NovaSeq X data poly-G tails (≥ 10)
  are on 0.06–0.31 % of reads and adapters on short-insert libraries much more
  ([certa/bench/results/hg002_2026-10/recent_illumina_costs.txt](certa/bench/results/hg002_2026-10/recent_illumina_costs.txt)).

### 3.6 Stage 5: output

- **Alignments**: SAM as now, or BAM, compressed by fastp-gpu's BGZF
  compressor pool (libdeflate) or on the GPU (nvCOMP deflate,
  `cuda_gzip_compress.cu`).
- **Fallback reads**: the uncertified reads as FASTQ.gz (the only compressed
  FASTQ left), or streamed uncompressed into a fallback aligner's stdin
  (`--fallback-cmd "bwa-mem2 mem ..."`), merged into the same BAM.
- **Failed reads**: fastp's `--failed_out`.
- **Reports**: fastp's JSON/HTML report and `certa`'s `--stats` JSON, plus a
  combined JSON with both.

### 3.7 Threads, streams and memory

- One CUDA context, created on the main thread before any worker starts
  (fastp-gpu initialises CUDA lazily because of a glibc/CUDA heap race seen
  during adapter detection; creating the context first, before any thread,
  avoids the race as well).
- Host threads: one reader per input file, a pool of workers that take whole
  batches in chunks (fastp's per-pack work, then CERTA's CPU passes on the same
  chunk, which `fused_chunk` already does for CERTA), and one writer per
  output with a compressor pool.
- GPU streams: input (inflate/parse), QC, mapping; consecutive batches
  overlap (batch n maps while batch n + 1 is QC'd).
- Device memory on an A100 80 GB: index and reference (GRCh38 `-s 8`:
  about 9 GiB; `-s 1`: 21.3 GiB), fastp's slots (1,360 MiB peak today), and
  two or three batches of 1 M reads (text about 0.4 GB each, codes 0.15 GB,
  results 0.2 GB).

## 4. Semantics to settle

1. **Trimming in the output** (soft clips by default; `--hard-trim`).
2. **Filtered reads** (low quality, too short, too many N): not mapped; to
   `--failed_out` and counted in both reports.
3. **Deduplication** (`--dedup`): duplicates are dropped before mapping, which
   saves mapping time (on 20 M HG002 NovaSeq X reads, 1.9 % of single reads
   and 0.24 % of pairs are exact duplicates).
4. **Merged pairs** (`-m`): mapped as single reads of up to 320 bases (CERTA's
   `LMAX`), with the certificate of a single read.
5. **Thread-dependent fastp statistics**: fastp-gpu's documented differences
   from upstream with more than one thread (insert-size histogram sampling,
   per-sequence adapter counts) stay as they are; nothing CERTA uses depends on
   them, except the `--bio` fragment-length model, which is built from a fixed
   prefix of the input instead.

## 5. Code organisation

A new repository, `certa-qc`, with both projects as submodules (both MIT,
fastp's upstream MIT too) and a CMake build:

- **fastp → `libfastp_core.a`** (built in `d0bromir/fastp`, published to
  fastp-gpu; M0). Everything except `main.cpp`: option parsing, the evaluator,
  the per-pack cores, statistics, reports and GPU kernels. The reader and
  writer threads stay in fastp's processors. Programs include one header,
  `fastp_core.h`, with its own namespace and plain types (fastp's internal
  headers put `Read`, `Options` and `Stats` in the global namespace).
  A `fastpcore::Session` takes fastp's command line, hands out batches,
  processes them on any worker slot, and writes fastp's reports. Each read
  comes back as a record with its front offset and modified flag.
- **CERTA**: `BatchSource` instead of `FastqInput` in `MapPipeline`,
  `GpuMapper::map_device`, soft clips for trimmed bases in `sam_output.cpp`.
- **`certa-qc`**: the driver (options of both programs, one `--threads`), the
  batch type and the bridge kernel.

```
certa-qc ref.cidx R1.fq.gz R2.fq.gz --gpu --device 1 -t 64 \
    [fastp: --detect_adapter_for_pe --trim_poly_g -q 15 -l 50 --dedup ...] \
    [certa: --intervals | --full-cert | --bio ...] \
    -o out.bam -u fallback_R1.fq.gz -U fallback_R2.fq.gz \
    --qc-json qc.json --qc-html qc.html --stats map.json
```

## 6. Validation

- **Equivalence.** For every option set of fastp-gpu's synthetic equivalence
  suite and certa's own test matrix: `certa-qc` output equals
  `fastp-gpu -o trimmed.fq` (uncompressed) followed by `certa map trimmed.fq`
  with `--hard-trim` (SAM byte for byte, fallback FASTQ, every report field
  fastp-gpu's suite compares), on 1, 8 and 64 threads, several batch sizes, CPU
  and GPU. Without `--hard-trim` the same alignments with the clips added.
- **Speed and memory.** Wall and CPU time and peak host/device memory against
  the two-step pipeline (fastp-gpu writing `.fq.gz`, then `certa map`) on
  fastp-gpu's eight benchmark datasets and 20 M HG002 NovaSeq X pairs, on
  galaxy (A100, ARM) and a2 (x86-64, CPU only), one session per host.
- **Placement accuracy.** Unchanged by construction (same reads, same mapper);
  checked through the equivalence test.

## 7. Milestones

| | Deliverable | Done when |
|---|---|---|
| M0 | `libfastp_core` split out of fastp (done on CPU, 2026-11-07) | fastp's extended equivalence suite matches the build before the split; the library's reads and reports match fastp's; GPU build still to test on galaxy |
| M1 | CPU `certa-qc`: ISA-L input, fastp per-pack work and CERTA's CPU passes in the same worker chunks, SAM output (built 2026-11-12 with two thread pools and hard trim; output equal, speed on 4 cores measured; see certa/README.md) | equivalent to the two-step pipeline; no intermediate file |
| M2 | GPU: one context, device batches, the bridge kernel, `GpuMapper::map_device` | equivalent; each batch's text uploaded once |
| M3 | nvCOMP / GDS input straight into device batches; BAM output (BGZF pool or nvCOMP) | equivalent; host never holds the text with GDS |
| M4 | combined report; data-driven fragment-length model for `--bio`; fallback aligner on a pipe | `--bio` pairs use the overlap histogram; one BAM with fallback reads |

## 8. Risks

- **fastp's structure.** Its processors own their reader and writer threads and
  read global options. M0 kept the threads in the processors and moved only the
  per-pack work. fastp's output is unchanged on CPU; the GPU build is not yet
  tested.
- **CPU contention.** fastp's per-pack work and CERTA's CPU passes share the
  cores, and `--full-cert` and `--bio` move most of CERTA's work to the CPU.
  The worker pool takes both per chunk, so neither starves, but the stage
  balance has to be measured (M1). First measurement, with two pools on a
  4-core machine (certa/README.md, "Speed on four cores"): certa map took
  80 % of the CPU time. With four threads each, `certa-qc` was 1–4 %
  slower than fastp, a plain file, then certa map, and 8–10 % faster than
  with a gzip file. With one thread each it was 8–12 % faster. So on the
  CPU the fused program gains where cores are free, not where mapping
  already fills them.
- **GPU memory.** With a `-s 1` index (21.3 GiB) and three batches in flight,
  an A100 40 GB is tight; batch size has to adapt to the free memory.
