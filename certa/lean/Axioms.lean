import Certa

/-
The headline theorems and the axioms they rest on (run by CI:
`lake env lean Axioms.lean`). A `sorry` only warns in `lake build`; here it
shows as `sorryAx`, and CI fails on it or on any axiom besides Lean's
standard three (propext, Classical.choice, Quot.sound).
-/
#print axioms Certa.lemmaL
#print axioms Certa.lemmaL_length_two
#print axioms Certa.span_eq
#print axioms Certa.lemmaE
#print axioms Certa.lemmaE_attained
#print axioms Certa.lemmaH
#print axioms Certa.lemmaH_length_three
#print axioms Certa.SpanS.lemmaS
#print axioms Certa.SpanS.lemmaS_length_four
#print axioms Certa.lemmaW
#print axioms Certa.diag_drift
#print axioms Certa.diag_within_pad
#print axioms Certa.no_split_within
#print axioms Certa.Greedy.greedy_max
#print axioms Certa.Greedy.greedy_max_tau
#print axioms Certa.Prune.prefix_live
#print axioms Certa.Prune.remove_pair
#print axioms Certa.Prune.above_live
#print axioms Certa.BandRow.vec_eq_scalar
#print axioms Certa.AffineFill.vec_eq_scalar
