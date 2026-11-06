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
     with other options (e.g. `--max-dist 1200`), called;
     [lost_reads.py](lost_reads.py) counts the reads minibwa places with
     MAPQ ≥ 20 that CERTA puts elsewhere or leaves unmapped.
