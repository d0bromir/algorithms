# Adapt CERTA's search, never its certified answer

CERTA can become adaptive without losing scientific validity if each adaptive decision is classified by what it touches. Some levers change only how work is executed: thread and batch scheduling, the CPU/GPU split, reordering reads for cache locality, memoizing duplicate reads with the name-hash tie-break replayed, and caching exact seed hit lists. These keep output byte-identical by construction, so any metric may drive them, even wall-clock timing. They are validated by a determinism matrix that also fuzzes the controller's decisions. Other levers change how hard CERTA searches a given read: pass budgets, and skipping or reordering passes for reads predicted to be uncertifiable. These cannot make a certificate wrong, but they change which reads are certified and fields such as MAPQ and the tier tag. They should be per-read rules frozen after calibration, so that a read's certified status depends on that read alone and stays invariant to threads, batches, devices, input order and sharding. Parameters that define the certified problem, such as the scoring scheme and the fixed 850-bp proper-pair distance, must never adapt online. The prior art shows the stakes. BWA-MEM's batch-estimated insert size made output depend on thread count until consortia pinned `-K`. Exhaustive mappers (GEM, RazerS, search schemes) show that adapting the search per read leaves guarantees intact. The potential gains are lopsided. Exact duplicates cap memoization at about 2% of single-end alignments, unless duplicate reads are disproportionately expensive. Pass V alone, by contrast, takes about 56% of a paired `--full-cert` run, and the 2.4–2.6% of reads that end uncertified may consume much of it. Two things are missing. The deciding measurement is per-read time per pass, split by final outcome. The deciding validity steps are freezing the controller before any held-out run, reporting four modes (static, adaptive, calibrated-static, oracle) with nonparametric confidence intervals from at least six runs, and adding input-order, sharding and controller-fuzzing checks to the existing byte-identity tests.

## Cross-read statistics broke BWA-MEM's thread invariance; per-read adaptation never did

### Batch-estimated state is the classic determinism failure

BWA-MEM infers pair orientation and the insert-size distribution "from a batch of reads". The default batch is **10,000,000 bases × the number of threads**, so changing `-t` changes the batches, the estimates and the output. The man page's only reproducibility advice is to fix the batch size with `-K` ([bwa.1 man page](https://github.com/lh3/bwa/blob/master/bwa.1)). In a 2017 issue, a user got different SAM output at `-t 2/6/10/16` but identical output when re-running at the same thread count. Heng Li's answer was to pass a large `-K` ([lh3/bwa issue #121](https://github.com/lh3/bwa/issues/121)). The fix spread through the ecosystem. The CCDG functional-equivalence standard requires the hidden option `-K 100000000` "to achieve deterministic alignment results" ([CCDG PipelineStandard](https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md)). Parabricks' GPU fq2bam needs the same `-K` as the CPU run to avoid small paired-end mismatches ([Parabricks fq2bam](https://docs.nvidia.com/clara/parabricks/latest/Documentation/ToolDocs/man_fq2bam.html)). bwa-mem3 states that "any equivalence claim on this page is implicitly a claim at a fixed batch size" ([bwa-mem3 equivalence](https://bwa-mem3.readthedocs.io/en/latest/whats-different/equivalence.html)).

The same lineage holds two further lessons. First, BWA-MEM breaks ties among repetitive hits with a hash of the read's position in the input (`hash_64(id+i)`), so reordering a FASTQ changes placements even at a fixed `-K` ([bwamem.c](https://raw.githubusercontent.com/lh3/bwa/master/bwamem.c)). Second, bwa-mem2 advises `-K 1000000` with its ERT index "to see the speedups" ([bwa-mem2 README](https://github.com/bwa-mem2/bwa-mem2)). That advice exposes the real design error: one knob sets both the scheduling batch (throughput) and the statistics epoch (science), so speed and reproducibility trade against each other. bwa-mem3 also documents a phantom `XS` tag that differs between AVX2 and AVX-512BW/NEON builds because query padding was not zeroed ([bwa-mem3 equivalence](https://bwa-mem3.readthedocs.io/en/latest/whats-different/equivalence.html)). Device identity can fail with no adaptivity involved at all.

Other tools avoided the trap in three ways:

- **No cross-read state.** Bowtie 2 reseeds its random number generator for every read from the read name, sequence, quality string and `--seed`, and its pairing bounds `-I/-X` are fixed parameters ([Bowtie2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)). minimap2's pairing code estimates no insert-size distribution, and it breaks ties with a hash of the read name ([minimap2 pe.c](https://github.com/lh3/minimap2/blob/master/pe.c); [minimap2.1](https://github.com/lh3/minimap2/blob/master/minimap2.1)).
- **Statistics reset at input-defined boundaries.** Strobealign estimates read length once, from the first 500 reads, and uses it to pick one of seven pre-optimized seed parameter sets (canonical lengths 50–400) ([strobealign README](https://github.com/ksahlin/strobealign)). The current code builds a fresh insert-size estimator for each chunk of 1,000,000 nucleotides, reseeds its random number generator per chunk and restores chunk order ([strobealign main.rs](https://github.com/ksahlin/strobealign/blob/main/src/main.rs)). My reading of the code is that output is thread-invariant but depends on chunk size; the project does not document this. Even the documentation is inconsistent: the paper calls the read-length estimate a median, while the README calls it an average and the code computes a mean ([Sahlin 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9753264)).
- **Sample once, then freeze.** DRAGEN runs up to 100,000 pairs through the aligner at the start and uses the resulting insert-length model for all remaining pairs ([DRAGEN v3.8](https://support-docs.illumina.com/SW/DRAGEN_v38/Content/SW/DRAGEN/GPipelineMeanIns_fDG.htm)).

BBMap reportedly keeps a running insert-size average per thread and is nondeterministic in multithreaded paired mode. That claim rests on a forum post only ([SEQanswers](https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/35652-bbmap-aligner-for-dna-rnaseq-is-now-open-source-and-available-for-download?p=290145)).

| Pattern of cross-read state | Example | Invariant to threads, batch size, device | Invariant to input order and sharding | Parameter it smuggles in |
|---|---|---|---|---|
| None | Bowtie 2, minimap2 pairing | Yes | Yes, if ties hash read content or name | None |
| Sample, then freeze | DRAGEN (100,000 pairs), strobealign (500 reads) | Yes | No: depends on the input prefix | Prefix size; needs an override flag such as `-r` |
| Epochs defined by the input | strobealign chunks, BWA-MEM at fixed `-K` | Yes | No | Epoch size (CCDG pins `-K`) |
| Epochs scaled by thread count | BWA-MEM default | No | No | Thread count |
| Running state per thread | BBMap (forum report) | No | No | The scheduler |

CERTA today is in the top row. Its proper-pair objective uses a fixed distance (`max_dist = 850` in `include/certa/pair.h`), not an estimated insert-size distribution. It breaks ties with an FNV-1a hash of the read name (`src/mapper.h`). Like Bowtie 2's and minimap2's tie-breaks, that hash depends on read content, not on input position as BWA-MEM's does. Every cross-read adaptation moves CERTA down the table. Its byte-identity claim survives only the first three rows, and the third only if the adaptation epoch is kept separate from the scheduling batch.

### Exhaustive mappers already adapt the search without weakening the guarantee

The exact-search literature keeps what is guaranteed separate from how it is found:

- **GEM** builds each filtering region one base at a time in the FM-index until the region's candidate count falls below a budget *t*. The pigeonhole principle then keeps the search exhaustive within the error bound, whatever partition results. Its fast modes skip reads that yield too few regions and align the rest exhaustively, mapping **about 0.2% fewer reads at twice the speed** ([GEM, Nature Methods 2012](https://doi.org/10.1038/nmeth.2221)).
- **RazerS 3** chooses its filter parameters from a minimum sensitivity stated by the user ([razers3 man page](https://dyn.manpages.debian.org/trixie/seqan-apps/razers3.1)).
- **Yara** offers `full`, `high` and `low` sensitivity settings ([yara_mapper man page](https://dyn.manpages.debian.org/buster/seqan-apps/yara_mapper.1)).
- **The Optimal Seed Solver** uses dynamic programming to pick, for each read, the least frequent set of non-overlapping seeds, cutting average seed frequency **3-fold** ([OSS](https://doi.org/10.1093/bioinformatics/btv670)).
- **Optimum search schemes** run **about 35× faster** than backtracking for two errors on 101-bp reads ([Kianfar et al.](https://arxiv.org/abs/1711.02035)).
- **Dynamic partitioning** adapts each pattern's partition at run time, for any search scheme, and still retrieves every occurrence ([Renders et al. 2021](https://doi.org/10.1016/j.isci.2021.102687)).

Learned indexes are a second output-preserving pattern. BWA-MEME raises seeding throughput by up to 3.45× and LISA by up to 4.5×, and both report output identical to BWA-MEM2 ([BWA-MEME](https://www.biorxiv.org/content/10.1101/2021.09.01.457579v1.full.pdf); [bwa-mem2 README](https://github.com/bwa-mem2/bwa-mem2)). The counterexample is Qtip. It trains a model on simulated reads that mimic the input sample and predicts MAPQ from it, so one field of every record depends on the whole sample ([Qtip](https://pmc.ncbi.nlm.nih.gov/articles/PMC5557537)).

CERTA's pass V already belongs to the first family. For reads still uncertified after pass 3, it chooses parts for each read: the rarest disjoint q-mers, optionally searched with one edit (project README). None of the exact tools surveyed adjusts search effort from statistics of *earlier* reads. That gap is CERTA's opportunity, and it is also where BWA-MEM's failure lives. GEM's fast mode is the right model for the lossy half: skipping hard reads costs coverage, not correctness, and the coverage lost must be reported as a number.

## FastQC's run-level metrics cannot route reads; CERTA's own pass statistics can

FastQC is a per-sample diagnostic. Its most "adaptive-looking" modules are computed on a sample of the run. The duplication and overrepresented-sequence modules track only sequences first seen in the **first 100,000 reads**, and reads longer than 75 bp are truncated to 50 bp ([FastQC duplicate sequences](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html)). Per-sequence quality warns when the most common mean read quality is below 27 ([FastQC per-sequence quality](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/3%20Per%20Sequence%20Quality%20Scores.html)). Adapter content warns above 5% and fails above 10% of reads ([FastQC adapter content](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html)).

Such aggregates can set run-level priors. fastp, for example, turns on polyG trimming from the machine ID, because on two-colour instruments "G means no signal" ([fastp README](https://github.com/OpenGene/fastp)). They say nothing about which read will be expensive. The features that work per read are the ones existing filters already use:

- N count: fastp drops reads with more than 5 N, and Bowtie 2's `--n-ceil` filters reads with more than 0.15 × length.
- Share of low-quality bases: fastp drops reads with more than 40% of bases below Q15.
- Read length.
- Complexity: fastp's score is the share of adjacent bases that differ, with a 30% threshold.
- A 3′ adapter hit.

([fastp README](https://github.com/OpenGene/fastp); [Bowtie2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml))

Low complexity affects correctness as well as cost. Li identified erroneous realignment in low-complexity regions as one of the two major sources of false variant calls ([Li 2014](https://arxiv.org/pdf/1404.0929)).

Most evidence that QC predicts mapping outcome concerns mapping rate on old data. Trimming raised the mapped share of 33-bp GA-era human RNA-seq from **72.2% to as much as 97.0%** ([Del Fabbro et al. 2013](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3871669/fullTextXML)). On modern data the picture is different. Soft-clipping aligners rescue adapters and low-quality bases, and untrimmed gene-level quantification is as accurate while the whole analysis runs up to an order of magnitude faster ([Liao & Shi 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7671312)). WGS germline calling with DeepVariant showed no difference between raw and trimmed reads across seven GIAB samples ([Barbitoff & Predeus 2024](https://f1000research.com/articles/13-506)). Aggressive trimming without a minimum length changed the expression estimates of more than 10% of genes through spurious mapping of short reads ([Williams et al. 2016](https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/s12859-016-0956-2)). For CERTA, that means trimming pushes reads toward the repeat tier rather than out of trouble. I found no study relating FastQC metrics to mapping *time*.

Prevalence data put a ceiling on routing by QC class. Vendor QC of NovaSeq X Plus NA12878 WGS reports:

- **0.00–0.01%** of pairs with more than 10% N bases
- **0.00%** of pairs with more than half their bases low-quality
- **0.12–0.17%** of pairs adapter-related

([Novogene](https://amea.novogene.com/about/news/reaching-maximum-465gb-lane-with-an-average-q30-95-novogenes-novaseq-x-plus-data-at-a-glance)). Illumina reports at least 85% of NovaSeq X bases in the Q40 bin, and more than 90% for a 10B HG002 PCR-free run ([Illumina Q40](https://www.illumina.com/science/genomics-research/articles/data-quality-q-scores.html)).

At those thresholds, the "bad read" classes are an order of magnitude rarer than the **2.39% of single-end reads** (100 − 97.61) that CERTA leaves uncertified even with `--full-cert`. So either most uncertified reads look ordinary and are really repeats, or they fall into QC classes only at much looser thresholds: a single N, a short adapter overlap, a polyG tail.

One provenance effect deserves a direct check. By default, BCL Convert trims adapters, pads reads trimmed below 35 bp with N up to 35 bp, and turns reads shorter than 22 bp into all-N reads ([BCL Convert settings](https://support-docs.illumina.com/SW/BCL_Convert/Content/SW/BCLConvert/SampleSheets_swBCL.htm); [Illumina KB 000009157](https://knowledge.illumina.com/software/bcl-convert/software-bcl-convert-reference_material-list/000009157)). In FASTQs produced this way, "short insert", "adapter" and "trailing Ns" describe a single population. Three signs in the HG002 set would confirm it: a spike at length 35, reads ending in N runs, and all-N reads. Those reads are candidates for cheap, provable early exits.

Published routing methods split cleanly into lossless and lossy ones:

- **Exact-match diversion (GenStore-EM).** Reads that match the reference exactly go to a cheap path. The authors argue that a strong hash makes this accuracy-neutral. The speedup over software baselines grows with the exact-match share, reaching **2.07–6.05×** on read sets that are 75–85% exact matches ([GenStore](https://arxiv.org/pdf/2202.10400)).
- **Seed-count gating (GenStore-NM).** Cheap chaining runs only on reads with few exact seeds, and reads with low chaining scores are dropped. This reaches **22–29×**, but only because 99.65% of that input did not align at all ([GenStore](https://arxiv.org/pdf/2202.10400)). The payoff tracks prevalence.
- **Pre-alignment filters.** Kernel speedups are large, end-to-end speedups small: Shouji cut the execution time of five aligners by up to 18.8×, and GateKeeper-GPU gave 1.4× end to end. SneakySnake treats a **0% false-reject rate as a requirement** ([Shouji](https://academic.oup.com/bioinformatics/article/35/21/4255/5421509); [GateKeeper-GPU](https://web3.arxiv.org/pdf/2103.14978); [SneakySnake](https://arxiv.org/pdf/1910.09020)).
- **Early rejection (GenPIP).** Low-quality nanopore reads are dropped early, as part of a 41.6× whole-system speedup "with negligible accuracy loss". The contribution of early rejection alone is not isolated ([GenPIP](https://arxiv.org/pdf/2209.08600)).

For CERTA, only filters with zero false rejects can preserve the certified set. Chaining-score or quality thresholds can reorder or route work. Used to skip passes, they shrink certified coverage.

The best predictors are ones CERTA already computes. Bowtie 2 treats a read as repetitive, and re-seeds it, when total seed hits divided by the number of seeds that aligned exceed 300 ([Bowtie2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)). GenStore-NM gates on the number of exact seeds. CERTA's first pass already returns comparable signals at no extra cost: part occurrence counts, how many parts hit, the best score so far, and how far that score is from a perfect one. Add four per-read features:

- N count
- number of bases outside the Q40 bin
- a hit on the 3′ adapter prefix
- a complexity score

Together these make a cheap feature vector computed per read, with no cross-read state. The Q40-bin count is informative because the bin's empirical error rate is about 10⁻⁴: a 150-base read made entirely of Q40-bin bases contains any error only about **1.5%** of the time (1 − (1 − 10⁻⁴)¹⁵⁰).

I found no published short-read method that predicts certifiability from read features. The closest are GenStore-NM's seed gate and a PeerJ preprint reporting F1 0.5–0.95 for run-length and complexity features predicting alignment quality. I saw that preprint only as a search snippet ([PeerJ Preprints 27428](https://peerj.com/preprints/27428)).

## Exact duplicates cap memoization near 2% unless they are the expensive reads

The project scanned its 20,001,073-read HG002 NovaSeq X set for exact duplicates:

- **1.94%** of single-end reads repeat an earlier read when the reverse complement counts as the same sequence.
- **1.46%** repeat an earlier read as written.
- **0.24%** of pairs repeat an earlier pair exactly.
- **3,580** distinct sequences occur 10 or more times.

If a hit cost nothing, the count-weighted ceiling would be 1/(1 − 0.0194) ≈ **1.020×** single-end and 1.0024× paired.

These rates are far below position-based duplication: a 9% mean in the 30× NovaSeq 6000 PCR-free 1000 Genomes set, and 13.7% ± 3.7% across GIAB samples, both measured with Picard ([Byrska-Bishop et al.](https://www.biorxiv.org/content/10.1101/2021.02.06.430068v1.full); [Baid et al.](https://www.biorxiv.org/content/10.1101/2020.12.11.422022v1.full)). That gap is expected. A single sequencing error breaks exact identity, while Picard only matches 5′ positions. fastp's duplicate definition requires every base to match for the same reason ([fastp README](https://github.com/OpenGene/fastp)).

Two cautions apply to these numbers. First, duplicate hits grow with the number of reads the memo sees. If the 20 M-read set is a subsample, a memo spanning a full 30× run would hit more often than these figures suggest. Second, the 3,580 high-multiplicity sequences account for at least 35,800 reads (0.18%). They are probably not ordinary reads but all-N reads, polyG tails, adapter dimers or satellite copies, and each may cost many times the median read. What matters is the time-weighted hit fraction, which comes from joining per-read time with multiplicity. No published value exists.

The prior art offers precedents, but not this design:

- Bowtie 2 keeps a cache of seed substrings that is "particularly effective when the input is sorted or otherwise grouped" ([bowtie2 aligner_cache.h](https://raw.githubusercontent.com/BenLangmead/bowtie2/master/aligner_cache.h)).
- CORA compresses identical k-mers across reads and maps each unit once against a reference homology table. Its acceleration "improves as the total size of the … read dataset increases" ([CORA](https://cb.csail.mit.edu/cora/); [MIT group page](https://cb.csail.mit.edu/gems/)).
- Masai searches the seeds of a whole batch jointly and is 2–4× faster than Bowtie 2 and BWA ([Masai](https://academic.oup.com/nar/article/41/7/e78/1068067)).
- Small-RNA pipelines collapse identical reads before mapping ([fastx_collapser](https://nf-co.re/modules/fastx_collapser)).

minimap2 and BWA-MEM map every read independently, with no result cache. Both break ties on something other than sequence: minimap2 on the read name, BWA-MEM on the read's input position ([minimap2 map.c](https://raw.githubusercontent.com/lh3/minimap2/master/map.c); [bwamem.c](https://raw.githubusercontent.com/lh3/bwa/master/bwamem.c)). No paper describes a whole-read memo keyed by sequence that keeps a best-mapper's output byte-identical.

CERTA can do it because each record is a function of the read's content plus a tie-break hash of the read name. The memo must store the complete result as it stands before the tie-break:

- every co-optimal locus, in canonical order
- the certificate bounds and the inputs to XQ
- the tier
- whatever the MAPQ computation uses

On a hit, the FNV-1a hash of the new read's name picks the reported locus exactly as a recomputation would. Caching only the winning locus would break identity.

The key must include every input the per-read computation reads: the sequence, its length, and the quality string if any filter or score uses it. Every hit must be verified against the full key. fastp tolerates a **0.01%** hash-collision rate ([fastp README](https://github.com/OpenGene/fastp)); in CERTA, any collision would be a correctness bug.

Sharing results between a read and its reverse complement accounts for the 0.48-point difference between 1.94% and 1.46%. It is safe only if the alignment computation is strand-symmetric, including traceback tie-breaks and locus order. Otherwise, key on the read as given. In paired mode, a per-mate memo of the single-end stages reaches 1.46–1.94% of mates, but the joint pair certificate still has to run for every pair.

The memo's most useful property is that a hit or a miss cannot change the output. Capacity, eviction policy, races between threads and batch boundaries therefore cannot affect it: any cache state yields the same bytes. That permits a bounded memo, which keeps peak memory at the level the project has just reached, the index plus a constant (21.3 GiB single-end, 21.5 GiB paired). An unbounded memo would not: even samblaster's positional signatures need about 20 MB per million pairs ([samblaster README](https://raw.githubusercontent.com/GregoryFaust/samblaster/master/README.md)). Optical and pad-hopping duplicates cluster on tiles: in one patterned-flow-cell library, 92% of duplicate pairs lay on a single tile ([QCFail](https://sequencing.qcfail.com/articles/illumina-patterned-flow-cells-generate-duplicated-sequences/)). Because FASTQs are written tile by tile, a window over recent reads may capture most hits. The reuse-distance distribution from the existing scan would size that window.

The bigger reuse opportunity is probably high-frequency seeds and tier R results rather than whole reads. Alu and L1 cover about 11% and 17% of the genome ([bioRxiv 2024](https://biorxiv.org/cgi/content/short/2024.09.10.612283v1)), and tier R takes 26% of CPU time in paired interval mode. A cache limited to high-occurrence parts would catch repeat-derived reads that never repeat exactly. No one has measured its hit rate in FASTQ order.

## Valid adaptive claims need a frozen controller and untouched data

Benchmarking guidance converges on a few rules that bite hardest on adaptive methods:

- Every tool gets the same tuning strategy, and tuning one's own method more than the competitors' is a source of bias ([Weber et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6584985/)).
- Papers introducing a method reported it winning **94.7%** of the time, against 64.3% in neutral comparisons ([Buchka et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8111726/)).
- Calibrating and validating on the same data is a benchmarking crime ([Heiser](https://gernot-heiser.org/benchmarking-crimes.html)).
- Data leakage has affected 294 papers in 17 fields ([Kapoor & Narayanan 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10499856/)).

Truth sets and their statistics are standard. GIAB v4.2.1 covers 92% of HG002's autosomal GRCh38 ([Wagner et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9706577/)). GA4GH practice compares calls with hap.py and vcfeval, stratifies by region and reports confidence intervals ([Krusche et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6699627/)). GIAB stratifications exist for GRCh38 and CHM13 ([Dwarshuis et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11489684/)). Holding out samples is community practice. Truth Challenge V2 published the HG002 truth and kept HG003 and HG004 blind ([NIST](https://www.nist.gov/node/1838346)), and DeepVariant v0.5 withheld HG002 and chromosome 20 from training ([Google Research](https://www.research.google/blog/deepvariant-accuracy-improvements-for-genetic-datatypes/)). Rabema's all-best category fits certified claims. Its authors also showed that a mapper's parameter settings alone can shift its sensitivity by more than 20 points ([Holtgrewe et al. 2011](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3128034/fullTextXML)).

For CERTA, the structural argument is strong, but it has to be split by lever. When a lever leaves output byte-identical, accuracy cannot overfit at all. Only speed can, so the held-out samples must show that the speed gains generalize. When a lever changes the certified set, the certified share and the content of uncertified records become tuned outcomes. The controller's thresholds, feature set, warm-up length and epoch size are then hyperparameters, subject to the hold-out rule.

The current protocol calibrates on HG002 and holds out HG001 and HG005, which mirrors Truth Challenge practice. If HG001 or HG005 informed the design of the existing passes, a learned router needs at least one sample that played no role in design. That could be another GIAB genome, ideally with a different instrument or read length, because GIAB samples share data-generation pipelines. The cleanest protocol freezes the controller as a commit hash before any held-out run. It then reports four modes for each sample and host:

- the static defaults
- the frozen adaptive controller
- the best static setting from HG002, applied unchanged
- an oracle static setting chosen on the test sample itself, labelled as an upper bound

Reviewers will look first for the gap between the adaptive mode and the calibrated-static mode on held-out data. Calibrating on one host and testing on the other adds hardware as a second hold-out axis.

Timing statistics need tightening too. Hoefler and Belli require absolute values next to every ratio, a statement of whether the measurements are deterministic, nonparametric confidence intervals, and no arithmetic means of ratios. In their survey, only 2 of 95 HPC papers reported confidence intervals ([Hoefler & Belli 2015](https://htor.inf.ethz.ch/publications/img/hoefler-scientific-benchmarking.pdf)). With three runs, the min–max range contains the true median only 75% of the time (1 − 2 × 0.5³), while six runs give a 96.9% nonparametric interval. Medians of three are defensible only if labelled as ranges, and headline adaptive-versus-static comparisons should use at least six interleaved runs.

Speed comparisons must also hold reproducibility constant. BWA-MEM2 should run at a fixed `-K`, as CCDG specifies. If CERTA adds a fast adaptive mode, the competitors should get their documented fast presets.

The determinism checks have to grow with the controller. The project already checks byte identity across builds, hosts and CPU vs GPU. An adaptive CERTA needs three more checks:

- **Input order:** shuffle the FASTQ and compare after sorting by read name.
- **Sharding:** split the input into *k* parts, map each separately and concatenate.
- **Controller fuzzing:** replace every execution decision (batch size, device, order) with draws from a seeded random number generator, over many seeds, and require an identical SHA-256 of the SAM body.

Certificate validity has its own arbiter. In certifying-algorithm practice, trust rests on an independent checker that checks the original problem, not just the witness's internal consistency ([McConnell et al.](https://researchportalplus.anu.edu.au/en/publications/certifying-algorithms/); [Beyer 2025](https://www.sosy-lab.org/research/prs/2025-06-03_Dagstuhl-CertAlg_Witnesses-2.0_Dirk.pdf)). SAT competitions work the same way, with DRAT proofs and a verified checker, cake_lpr ([DRAT-trim](https://arxiv.org/pdf/1610.06229); [cake_lpr](https://research.chalmers.se/publication/531575/file/531575_Fulltext.pdf)). Every lossy mode should pass `certa check` with zero rejections on every record of the held-out runs. That turns "no certificate became wrong" from an argument into a measurement.

## Three classes of levers need three different evaluations

The project's measurements set the stakes. `--full-cert` raises the certified share from **94.38% to 97.61%** single-end and from **95.39% to 97.44%** paired, at about **2.5× the default time**. The README gives 38.7 s single-end and 100.1 s paired on 64 cores. Within those runs:

- pass V is about **56%** of a paired `--full-cert` run
- banded dynamic programming is about **44%** of a single-end `--full-cert` run
- tier R is about **26%** of CPU time in paired interval mode

Passes 3, V and L run on the CPU. The extra 3.23 points single-end therefore cost about 1.5 default runs. Meanwhile 2.39% of single-end and 2.56% of paired reads end uncertified, after passing through every expensive stage without changing status.

Suppose a fraction X of pass V's time goes to reads that end uncertified. A perfect skip would then save 0.56·X of a paired `--full-cert` run. At X = 0.5 that is 28%, a 1.39× speedup. This is only an illustration, because X has not been measured. One measurement decides whether routing or memoization deserves the investment: per-read time per pass, split by final outcome and joined with duplicate multiplicity.

| Lever | Driven by | Effect on output | Condition that must hold | How to evaluate it |
|---|---|---|---|---|
| Thread count, work stealing, adaptive batch size | Timing, queue depth, memory | Byte-identical | No statistic computed per batch; output order restored; no per-read state carried in threads | Determinism matrix plus seeded fuzzing of batch sizes; SHA-256 of SAM body |
| CPU/GPU split, including sending predicted-hard reads to the CPU passes early | GPU occupancy, per-read cost prediction | Byte-identical | Every kernel bit-identical on every read class; integer scoring; zeroed padding | Fuzzed split; CPU-vs-GPU diff broken down by tier |
| Reordering reads for locality (e.g. minimizer buckets) | Seed statistics | Byte-identical | Ties broken by name hash, never by input position; reorder buffer on output | Shuffled input, diff after sorting by name |
| Duplicate memo with tie-break replay | Sequence multiplicity | Byte-identical | Full key including any quality used; full canonical tie set; full-key verification; reverse-complement sharing only if strand-symmetric | Memo off, on, and with one entry; hit rate, time-weighted saving, GiB |
| Cache of high-frequency seeds and tier R results | Seed occurrence | Byte-identical | Cached object equals the exact index or tier R result | As for the duplicate memo |
| Skip when futility is provable | The read's own bounds from earlier passes | Same certified set; bytes identical only if a failed pass leaves the record unchanged | A lemma that the pass cannot certify the read | Zero diff against a run without skipping, on all samples |
| Reordering passes per read | Read features | Certified status can change; tier, alignment and MAPQ can change | Fields computed canonically after certification | Semantic diff (POS, CIGAR, AS, XL = XU) plus byte diff |
| Skipping reads predicted uncertifiable (frozen per-read rule) | N, adapter, complexity, non-Q40 bases, seed counts, score gap | Certified set shrinks by the false skips | Rule frozen and versioned in `@PG`; audit sample selected by hash | Certified-share loss with CI, speedup, Pareto curve, loss by stratum, held-out samples |
| Budgets adapted from earlier reads over input-defined epochs | Yield of passes in earlier epochs | Certified set depends on earlier reads and on epoch size | Epochs counted in input reads; lag ≥ 1; integer reduction in fixed order; trace in header; freeze file | Everything above, plus shuffle and shard tests (expected to differ; quantify the difference) |
| Parameters that define the problem (scoring, `max_dist`, quality-aware penalties) | None online | Changes what is being certified | Set only in a profiling phase and declared | Evaluate as a different method |

### Execution levers can use any signal

The first five rows change only execution, so they may use non-deterministic signals such as timing and GPU queue depth, provided the output is proven, and tested, to be independent of every decision they make. This inverts BWA-MEM's mistake, where a throughput knob silently fed a scientific estimate.

Two of these levers carry a hidden risk. A dynamic CPU/GPU split, or a memo, turns any latent discrepancy into run-to-run nondeterminism. A read that lands on a different device, or is computed from a different copy, in different runs exposes any bit-level difference between the two paths, as bwa-mem3's SIMD-dependent `XS` tag did. The existing CPU-vs-GPU identity check is therefore a prerequisite. It should be broken down by tier, so that rare read classes are covered. Controller fuzzing is the new test.

Evaluation reduces to identity plus speed. The speed evidence should be medians of at least six runs with nonparametric intervals on both hosts, the overhead of each lever when it is enabled but idle, and speed gains generalizing from HG002 to the held-out samples.

### Search-effort levers need loss accounting

The middle rows change how hard CERTA searches, so they need loss accounting even when their certificates remain sound. Three details matter.

First, a "lossless" futility skip is byte-identical only if a failed pass never tightens XL/XU and never finds a better alternative locus. The README says default MAPQ is a point estimate that uses "the best other locus found", so the effort spent on a read can change its MAPQ even when its certified alignment does not.

Second, reordering passes is not neutral. A read sent straight to pass L would be certified as a local alignment (tier SL) where pass V might have certified an end-to-end one. `interval.h` notes that S0/S1 certify end-to-end edit distance, not the local score, so the two records can differ in their alignment, not only in their tag.

Third, a frozen per-read router keeps every invariance CERTA claims today, because a read's status depends only on the read and the rule. It gives up only equality with the static mode, and it can be made self-auditing. A deterministic audit subset, selected by a hash of the read name (for example 1% of reads), always runs the full pipeline. The stats file then reports the estimated certified share lost by skipping, with a binomial confidence interval, on every production run. That estimate is invariant to threads, input order and sharding.

### Cross-read adaptation and problem parameters

Adapting from earlier reads over input-defined epochs buys little beyond a per-read rule. It also costs input-order and shard invariance, and it makes a read's certified status depend on its neighbours in the FASTQ. That is the option to offer last, opt-in and with a trace.

For sharded pipelines, a better form of sample-then-freeze is a separate profiling step. Run it on a prefix or a hash-selected sample, have it write a parameter file whose hash goes in the SAM header, and pass that file to every shard, much as strobealign's `-r` makes its read-length estimate explicit. Parameters that define the problem stay fixed. Adapting them would change what CERTA proves, not how quickly it proves it.

## Conclusion

The adaptive question for an exact aligner differs from the one BWA-MEM faced. BWA-MEM needed cross-read statistics for accuracy, and the batch-level insert-size model was part of its answer. CERTA's certificates need no statistics from other reads, so cross-read state enters only for speed. The levers with the largest potential payoff can be placed either outside the output entirely (scheduling, device split, memo, seed caches) or inside a frozen per-read rule (routing on pass-1 statistics). Online adaptation of search effort from earlier reads is therefore mostly unnecessary. Where it is used, it introduces a new kind of bias that heuristic aligners never had to name. If certified status depends on a read's neighbours, two samples processed along different adaptive trajectories differ in which reads are certified. That difference falls in exactly the difficult regions where certificates matter most, and comparisons across samples would inherit it. Keeping certified status a function of the read alone avoids this problem by construction.

The work should proceed in a clear order. First, measure time per pass by final outcome and by duplicate multiplicity; this costs little and decides whether routing or memoization is worth building. Next, ship the byte-identical levers, with controller fuzzing added to the identity matrix. Then offer a frozen, hash-audited router as a named mode, with its certified-share loss reported on samples that played no role in its design. Each step adds a claim CERTA can test rather than assert: identical bytes for the first class, a measured and bounded coverage loss for the second, and zero checker rejections for both.
