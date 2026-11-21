# Variant calling with CERTA's MAPQ options

The scripts of [VALIDATION_RESULTS.md](../VALIDATION_RESULTS.md), "Variant
calling" (2026-11-05), as they ran. galaxy (ARM64) maps the reads, a2
(x86-64) calls and scores. The two hosts cannot reach each other, so a
third machine relays the records as CRAM.

1. **a2: tools and truth.**
   - [setup.sh](setup.sh) downloads the GIAB truth sets (HG002 v5.0q and
     CMRG v1.00, HG005 v4.2.1) and the v3.6 stratifications, and creates
     hap.py's environment.
     - Its bioconda `deepvariant=1.10.0` environment is not used: the
       package holds only wrapper scripts with unfilled placeholders.
   - DeepVariant comes from its Docker image instead. a2 runs no
     containers, so:
     - `python3 -I dv_fetch.py google/deepvariant 1.10.0 ~/vc/dv` unpacks
       the layers, checking each digest;
     - [dvrun.sh](dvrun.sh) runs each tool with the image's own loader and
       Python ([dv_unzip.py](dv_unzip.py) restores the symlinks of its zips).
   - [truth.sh](truth.sh) renames the truth sets' contigs to the
     reference's, writes hap.py's stratification list and RTG's SDF.
2. **galaxy: mapping.** [vc_map.sh](vc_map.sh) `GENOME OUT_DIR` streams the
   30× FASTQs once to five mappings ([fanout.py](fanout.py)). Each mapping
   keeps the records overlapping `targets.bed`: chr20 and the CMRG genes,
   padded by 1 kb.
3. **Relay.** [relay.sh](relay.sh) `GENOME`: CRAM on galaxy, a copy to a2,
   an md5 check, then step 4.
4. **a2: calling and scoring.** [call.sh](call.sh) `GENOME NAME...` runs
   DeepVariant (64 make_examples shards, call_variants,
   postprocess_variants) and hap.py with vcfeval. [vc_table.py](vc_table.py)
   `DIR GENOME EVAL NAME=LABEL...` writes the tables, with each run's
   mapping cost.
5. **Diagnostics.**
   - [remap60.sh](remap60.sh): every MAPQ ≥ 30 set to 60.
   - [vc_subset.sh](vc_subset.sh) `GENOME MODE...`: the same read pairs
     mapped again with each of CERTA's MAPQ modes.
   - [allele_mapq.py](allele_mapq.py): at SNV sites, the MAPQ of the reads
     that show the variant against those that show the reference;
     [allele_tags.py](allele_tags.py): their tiers.
   - [neutral.sh](neutral.sh), [neutral3.sh](neutral3.sh) (galaxy): `--mapq
     pair` models fitted without some inputs (`truth_mapq_pair.py fit
     --drop`), cross-validated on HG002; [vc_model.sh](vc_model.sh) (a2)
     builds CERTA with such a model and calls the subset.
   - [vc_opts.sh](vc_opts.sh) `GENOME TAG OPTIONS...`: the subset mapped
     with other options (e.g. `--max-dist 1200`), by `$CERTA` if set,
     called; [lost_reads.py](lost_reads.py) counts the reads minibwa
     places with MAPQ ≥ 20 that CERTA puts elsewhere or leaves unmapped,
     and why those with a nearby mate were not rescued.
   - [mate_rescue.sh](mate_rescue.sh) `BUILD_DIR`: `--iv-rescue local`
     against the same build without it, on both subsets
     ([VALIDATION_RESULTS.md](../VALIDATION_RESULTS.md#a-local-mate-rescue)).
     [swap_lost.py](swap_lost.py) gives the reads CERTA loses minibwa's
     records (with `--mapq`, also minibwa's MAPQ where CERTA places a read
     at its locus): how much of the gap to minibwa they explain.
   - [pr_mapq.sh](pr_mapq.sh) `BUILD_DIR` (a2): the subsets with PR mates'
     MAPQ from the pairs found; [same_records.py](same_records.py) checks
     that only MAPQ and XQ differ from the build before
     ([VALIDATION_RESULTS.md](../VALIDATION_RESULTS.md#the-mapq-of-mates-placed-by-the-pair)).
   - [t2t_feat.sh](t2t_feat.sh) `BUILD` (galaxy): the subsets mapped with
     `--t2t --feature-tags`, one TSV of T2T classes per genome (relayed to
     a2's `~/vc/diag`); [t2t_gate.py](t2t_gate.py) (a2) writes the BAMs of
     the T2T-gated MAPQ rules for call.sh. [errors_table.py](errors_table.py)
     `DIR EVAL NAME=LABEL...`: FN and FP per run, in all, segdups and the
     difficult regions.
   - [lost_run.sh](lost_run.sh) (a2): the reads minibwa places with MAPQ
     ≥ 20 that CERTA puts elsewhere, re-mapped with `--feature-tags
     --witness`; [lost_diag.py](lost_diag.py) re-scores minibwa's alignment
     under CERTA's objective against CERTA's interval and parts (search
     miss, scored as well, or above XU).
     [vote_sweep.sh](vote_sweep.sh): the lost pairs with seed votes off and
     at several settings ([VALIDATION_RESULTS.md](../VALIDATION_RESULTS.md#reads-certa-lost)).
   - [dv_map.sh](dv_map.sh) `BUILD_DIR TAG OPTIONS...` (galaxy): the
     subsets mapped there (a2's disk is nearly full) into `S_TAG.bam`;
     [dv_relay.sh](dv_relay.sh) `TAG` (the relay machine) copies them to a2
     and calls them.
   - [bwscale.sh](bwscale.sh) `GENOME IN OUT TABLE` (a2): a BAM or CRAM with
     every MAPQ rewritten through a table, for the fitted MAPQ on BWA-MEM's
     scale ([VALIDATION_RESULTS.md](../VALIDATION_RESULTS.md#the-fitted-mapq-on-bwas-scale)).
   - [dv_final.sh](dv_final.sh) `BUILD_DIR TAG=OPTIONS...` (galaxy): the
     subsets mapped with several option sets, the reads extracted once, as
     CRAM (a third of the BAM, for the 3 MB/s relay);
     [dv_relay_cram.sh](dv_relay_cram.sh) `TAG...` (the relay machine)
     copies them to a2, checks their md5 and starts call.sh there.
   - Error sites, read by read (a2; the sites as TSV from hap.py's VCFs):
     - [site_reads.py](site_reads.py): the reads over indel sites per run,
       with the indel or clipped near it
       ([hg005_indel_sites.md](../results/truth_2026-11_refit/hg005_indel_sites.md)).
     - [snv_tiers.py](snv_tiers.py), [snv_missing.py](snv_missing.py) and
       [snv_ties.py](snv_ties.py): at false SNV calls only CERTA makes, the
       reads with the false allele by tier and what minibwa does with them,
       the reads minibwa has there that CERTA lacks, and their equally best
       loci.
     - [lowend.sh](lowend.sh): the model's Phred 3–4 rewritten below
       DeepVariant's MAPQ floor and called
       ([snv_fp_diagnosis.md](../results/truth_2026-11_fp/snv_fp_diagnosis.md)).
     - [snv_t2t.py](snv_t2t.py): where T2T-CHM13 puts the reads CERTA counts
       at those sites (from a `--t2t --feature-tags` run of the same reads,
       [t2t_feat.sh](t2t_feat.sh)); [t2t_tiers.py](t2t_tiers.py): T2T's
       disagreements over all records, by tier and region
       ([lean_t2t.md](../results/truth_2026-11_lean/lean_t2t.md)).
