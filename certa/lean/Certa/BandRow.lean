/-
The vectorized DP row (THEORY.md "The band DP row"; certa/longread.h
band_row against band_row_scalar). With h0[b] = max(m, ins) of band cell b
(both computed from the previous row, the same in both versions), the
scalar recurrence

  d = max(left − 7, d − 1),  h[b] = max(h0[b], d),  left = h[b]

started from the sentinel NEG (left = d = NEG), and the vector form

  incl = max(carry, h0[b] + b),  h[b] = max(h0[b], incl − b − 6),  carry = incl

started from carry = NEG (the prefix maximum over b' ≤ b, computed by
shifts within a vector and the carry across vectors) give the same row
whenever every h0 is at least NEG: every h0 is a real score (≥ −9), NEG is
INT_MIN / 4. Here for rows of any width.
-/
namespace Certa.BandRow

def scalarAux : Int → Int → List Int → List Int
  | _, _, [] => []
  | left, d, x :: t =>
    let d' := max (left - 7) (d - 1)
    let h := max x d'
    h :: scalarAux h d' t

def vecAux : Int → Int → List Int → List Int
  | _, _, [] => []
  | carry, b, x :: t =>
    let incl := max carry (x + b)
    max x (incl - b - 6) :: vecAux incl (b + 1) t

/-- The invariant: the scalar deletion value about to be formed at b and the
vector carry agree up to the sentinel's own decay. -/
theorem rows_eq (NEG : Int) :
    ∀ (h0 : List Int) (left d carry b : Int),
      (∀ x ∈ h0, NEG ≤ x) → 0 ≤ b →
      max (left - 7) (d - 1) = max (NEG - 1 - b) (carry - b - 6) →
      scalarAux left d h0 = vecAux carry b h0
  | [], _, _, _, _, _, _, _ => rfl
  | x :: t, left, d, carry, b, hx, hb, hinv => by
    simp only [scalarAux, vecAux]
    have hxN : NEG ≤ x := hx x (by simp)
    have ht : ∀ y ∈ t, NEG ≤ y := fun y hy => hx y (by simp [hy])
    have heq : max x (max (left - 7) (d - 1)) = max x (max carry (x + b) - b - 6) := by
      rw [hinv]; omega
    have hrest := rows_eq NEG t (max x (max (left - 7) (d - 1))) (max (left - 7) (d - 1))
      (max carry (x + b)) (b + 1) ht (by omega) (by rw [hinv]; omega)
    rw [hrest, heq]

/-- **The vector row equals the scalar row** for any row of values above
the sentinel. -/
theorem vec_eq_scalar (NEG : Int) (h0 : List Int) (hx : ∀ x ∈ h0, NEG ≤ x) :
    scalarAux NEG NEG h0 = vecAux NEG 0 h0 :=
  rows_eq NEG h0 NEG NEG NEG 0 hx (Int.le_refl _) (by omega)

end Certa.BandRow
