# Machine-checked proofs (Lean 4)

Lean 4 proofs of the lemmas behind CERTA's bounds ([THEORY.md](../THEORY.md)).
`tests/verify_lemmas.cpp` checks the same statements exhaustively up to a size;
here they hold for every size. The toolchain is pinned in `lean-toolchain`
(Lean 4.35.0); no Mathlib.

    lake build    # checks every proof (about 2 s)

Install the toolchain with [elan](https://github.com/leanprover/elan); `lake`
then picks the pinned version.

| File | Statement | THEORY.md | Finite check |
|---|---|---|---|
| `Certa/LemmaL.lean` | `lemmaL`: every part of length ≥ 3 that does not occur exactly costs ≥ 5, disjointly, for transcripts of any length. `lemmaL_length_two`: parts of length 2 do not suffice. | §4, Lemma L | `verify_lemmas`: all transcripts of reads up to 9 bases |

The model is `verify_lemmas.cpp`'s: each read base is a match, a mismatch or
inserted, and may be followed by a deletion; the loss is 5 per mismatch, 2 per
inserted base plus 6 per run of them, and 6 + k per deletion of k. The proof
credits every base with its own loss, splitting a run's opening cost of 6 as 3
to its first base and 3 to its last; the credits sum to at most the loss, and
every inexact part of length ≥ 3 collects at least 5.
