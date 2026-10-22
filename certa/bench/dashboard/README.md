# CERTA vs competitors: the dashboard's source

The published dashboard (a claude.ai Dashboard artifact) is the page
`index.html` over the datasets in `data/`, every one derived from
`bench/results/`:

| dataset | from | made by |
|---|---|---|
| `speed.csv` | `results/hg002_2026-10/final_{galaxy,a2}.tsv` | `make_speed_csv.py results data/speed.csv` |
| `slow_speed.csv`, `slow_reliability.csv`, `slow_hifi.csv` | `results/hg002_2026-10/slow_*` | `make_slow_csv.py results data` |
| `provable_errors.csv` | `results/hg002_2026-09/weaknesses_2M.txt` | by hand |
| `certification.csv` | the stats JSON files in `results/hg002_2026-09` and `hg002_2026-10` | by hand |
| `variant_calling.csv` | `results/hg002_2026-09/giab_chr20_vcfeval.tsv`, `results/hg002_2026-10/giab_chr20_vcfeval_*.tsv` | by hand |
| `hifi_cost.csv`, `hifi_outcome.csv` | `results/aim6_2026-10/` | by hand |

The page draws every figure from these datasets (no number is written in its
markup); the ratios it shows are computed in the page (`dash.calc`). Run the
scripts from `certa/bench`.
