# certa-qc against fastp, then certa map (2026-11-12)

`bench/qc_speed.py` with its defaults, run three times:

```bash
python3 bench/qc_speed.py --certa build/certa --certa-qc build/certa-qc --fastp ~/fastp/fastp --workdir WORKDIR
```

- **Data.** A 2 Mbp simulated genome (3 contigs) and 400 000 pairs of 150 bases with NovaSeq X
  qualities, gzip input. `tests/qc_equivalence.py` adds adapter read-through, poly-G tails and poor
  ends, as in the equivalence test.
- **Machine.** A virtual machine with 4 vCPUs (Intel Xeon, 2.10 GHz) and 15 GiB of memory. The
  inputs and the intermediate files stay in the page cache.
- **Builds.** `certa` and `certa-qc` from commit cf6b460 (the sources are unchanged since), Release,
  `-march=native`, GCC 13.3. `fastp` and `libfastp_core.a` from d0bromir/fastp commit 5b7d127
  (branch `claude/m0-libfastp-core`).
- **Script.** `run1.txt` and `run2.txt` come from the script at commit 79d5f50. `run3.txt` comes
  from commit 7f73f03, which also prints the two-step pipeline's fastp and certa map times.

Every run wrote the same SAM. Medians of three rounds per run (wall time / CPU time, in seconds):

| Pipeline | fastp `-w` | certa `-t` | run 1 | run 2 | run 3 |
|---|---|---|---|---|---|
| fastp, gzip file, certa map | 4 | 4 | 10.64 / 37.23 | 10.63 / 35.86 | 10.99 / 37.37 |
| fastp, plain file, certa map | 4 | 4 | 9.50 / 33.24 | 9.62 / 32.01 | 9.75 / 32.27 |
| certa-qc | 4 | 4 | 9.59 / 33.33 | 9.80 / 34.18 | 10.16 / 35.19 |
| certa-qc | 1 | 4 | 10.43 / 34.46 | 10.21 / 33.53 | 11.15 / 37.10 |
| certa-qc | 1 | 3 | 11.29 / 33.77 | 11.23 / 33.39 | 12.14 / 36.42 |
| fastp, plain file, certa map | 1 | 1 | 31.12 / 33.07 | 31.07 / 32.89 | 30.72 / 32.56 |
| certa-qc | 1 | 1 | 27.46 / 34.27 | 28.53 / 35.00 | 27.81 / 34.76 |
