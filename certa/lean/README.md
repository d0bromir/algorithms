# Machine-checked proofs (Lean 4)

Lean 4 proofs of the lemmas behind CERTA's bounds ([THEORY.md](../THEORY.md)).
`tests/verify_lemmas.cpp` checks the same statements exhaustively up to a size;
here they hold for every size. The toolchain is pinned in `lean-toolchain`
(Lean 4.35.0); no Mathlib.

    lake build                    # checks every proof (about 2 s)
    lake env lean Axioms.lean     # the axioms of every headline theorem

Install the toolchain with [elan](https://github.com/leanprover/elan); `lake`
then picks the pinned version. CI (the `lean` job of
`.github/workflows/certa.yml`) runs both, and fails on a `sorry` (which
`lake build` only warns about) or on any axiom besides Lean's standard three
(`propext`, `Classical.choice`, `Quot.sound`).

| File | Statement | THEORY.md | Finite check |
|---|---|---|---|
| `Certa/LemmaL.lean` | `lemmaL`: every part of length ≥ 3 that does not occur exactly costs ≥ 5, disjointly, for transcripts of any length. `lemmaL_length_two`: parts of length 2 do not suffice. | §4, Lemma L | `verify_lemmas`: all transcripts of reads up to 9 bases |
| `Certa/LemmaH.lean` | `lemmaH`: with parts of length ≥ 4, every part no one-mismatch probe witnesses costs ≥ 7, disjointly, for transcripts of any length (`Witnessed` as `verify_lemmas`' `witnessed1`). `lemmaH_length_three`: parts of length 3 do not suffice. | §11 list, Lemma H | `verify_lemmas`: every layout of 10-base reads; three parts of 4 at L = 12 |
| `Certa/LemmaS.lean` | `lemmaS`: for disjoint parts of length ≥ 5, Lemma L's bound over a span [x, y) is at most the full span's, for any counted set and extra (the second claim of Lemma S; it makes mutant L6 equivalent). `lemmaS_length_four`: parts of length 4 do not suffice. | §4, Lemma S | `verify_lemmas`: every layout of 24-base reads, every counted set, extra ≤ 2 |
| `Certa/TheoremG.lean` | `greedy_max`: for a family of read intervals closed under extension, the greedy partition (the shortest interval at each start, then on from its end) is a maximum set of disjoint intervals of the family. `greedy_max_tau`: the same for the intervals occurring ≤ τ times in a reference text, whose closure under extension is proved (`occ_super`). | §6, Theorem G | `verify_lemmas`: greedy = DP optimum on 20,000 random instances |
| `Certa/Prune.lean` | The long-read band pruning (`band_best` with a floor), for alignment paths: `prefix_live`, every cell on a path to a live cell is live; `above_live`, once no new piece can lead above the floor, the cell above a live cell is live (`remove_pair`: dropping an aligned pair after the first keeps the diagonal and loses at most 1). | "Pruning the band" | `test_longread`: the pruned DP equals the full one above the floor |
| `Certa/LemmaW.lean` | `lemmaW`: an alignment scoring above b has gap total ≤ max(0, L − b − 7), the pad of every certificate's band. | §4, Lemma W | — |
| `Certa/LemmaC.lean` | Lemma C's geometric core: `diag_drift`, two aligned bases lie on diagonals at most the gap total apart; `diag_within_pad`, at most pad(b) apart above b; `no_split_within`, no chain split falls between two diagonals at most pad apart. | §4, Lemma C | — |
| `Certa/BandRow.lean` | `vec_eq_scalar`: the vectorized DP row (an inclusive prefix maximum with a carry) equals the scalar Gotoh recurrence, for rows of any width, when every value is above the sentinel. | "The band DP row" | `test_longread`: exact equality on 20,000 random rows |
| `Certa/LemmaE.lean` | `lemmaE`: a local alignment whose completion has d ≥ 1 edits scores ≤ edit_bound(L, d) = L − d − (d ≥ 2 ? 5 : 4), for reads of any length; the score is BWA-MEM's (+1, −4, gaps −(6 + k), −5 per clipped end), equal to span − loss (`span_eq`). `lemmaE_attained`: the bound is attained for every d ≥ 1. | §4, Lemma E | `verify_lemmas`: all clipped transcripts of 8-base reads |

The model is `verify_lemmas.cpp`'s: each read base is a match, a mismatch or
inserted, and may be followed by a deletion; the loss is 5 per mismatch, 2 per
inserted base plus 6 per run of them, and 6 + k per deletion of k. Lemma L's
proof credits every base with its own loss, splitting a run's opening cost of 6
as 3 to its first base and 3 to its last; the credits sum to at most the loss,
and every inexact part of length ≥ 3 collects at least 5. The same credits give
every unwitnessed part of length ≥ 4 at least 7 (Lemma H). Lemma E's proof
carries the loss beyond the edits through the transcript.
