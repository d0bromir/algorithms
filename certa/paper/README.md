# The CERTA paper

A draft of an Original Paper for *Bioinformatics* (OUP). It uses the
article class, laid out like the journal's (structured abstract,
Introduction / Methods / Results / Discussion, author–year citations). The
journal's own class is applied at acceptance.

## Every number comes from the data

The text quotes no hand-typed measurement. Each number is a LaTeX macro in
`numbers.tex`, which `make_numbers.py` writes from the raw results in
`../bench/results/`, the mutant list `../tests/mutants.tsv` and the Lean
sources `../lean/`. Each macro's line in `numbers.tex` names its source
file. The figures are drawn by `make_figures.py` from the same files.

## Build

```bash
python make_numbers.py && python make_figures.py && latexmk -pdf main.tex
```

- **Needs:** Python 3 with matplotlib, and a LaTeX distribution with
  amsmath, mathptmx, natbib, booktabs, multirow, caption, enumitem, lineno,
  hyperref, authblk and TikZ.
- **With MiKTeX:** pass `-pdflatex="pdflatex --disable-installer %O %S"` to
  latexmk, so a missing package fails the build instead of opening an
  install prompt.
- **What is committed:** `numbers.tex` and `figures/*.pdf`, so the paper
  also builds with LaTeX alone. Re-run both scripts after any benchmark
  changes.

## Files

| File | What it is |
|---|---|
| `main.tex` | the paper |
| `refs.bib` | references, each checked against its Crossref, DataCite or arXiv record |
| `make_numbers.py` | results → `numbers.tex` |
| `make_figures.py` | results → `figures/` (palette: categorical slots validated for colour-vision deficiency) |

## Before submission

- Add the affiliation, the contact address and funding.
- Decide what goes into the supplement. `../THEORY.md` holds the written
  proofs and `../FINDINGS.md` the full measurements; both are cited as
  online material.
- Bring the length to the journal's limit (about 7 pages in its format).
