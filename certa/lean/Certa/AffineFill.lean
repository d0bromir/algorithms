/-
The vectorized traceback DP (THEORY.md "The traceback DP, vectorized";
certa/core.h affine_fill_vec against affine_fill_scalar). Within one row,
with h0[b] = max(m, E) of band cell b (computed from the previous row, the
same in both versions), the scalar loop forms the deletion state F, the
cell value h and the traceback bit "F extends F" by

  b = 0:  F = NEG,  h = max(h0, NEG),  bit = false
  b ≥ 1:  F = max(h[b-1] - 7, F[b-1] - 1),  h = max(h0, F),
          bit = F[b-1] - 1 > h[b-1] - 7

and the vector form by an exclusive prefix maximum with a carry that starts
at NEG + 6:

  b = 0:  F = NEG,  h = max(h0, NEG),  bit = false,  carry = max(NEG + 6, h0)
  b ≥ 1:  F = carry - b - 6,  h = max(h0, F),  bit = F[b-1] + 6 > h[b-1],
          carry = max(carry, h0 + b)

The two give the same F, h and bit in every cell, for any integers: unlike
the band DP row (BandRow.lean) there is no premise about the sentinel.
-/
namespace Certa.AffineFill

/-- The scalar loop from b = 1 on, with F[b-1] = f and h[b-1] = hl. -/
def scalarRest : Int → Int → List Int → List (Int × Int × Bool)
  | _, _, [] => []
  | f, hl, x :: t =>
    let F := max (hl - 7) (f - 1)
    let h := max x F
    (F, h, decide (f - 1 > hl - 7)) :: scalarRest F h t

def scalarRow (NEG : Int) : List Int → List (Int × Int × Bool)
  | [] => []
  | x :: t => (NEG, max x NEG, false) :: scalarRest NEG (max x NEG) t

/-- The vector form from b = 1 on: the carry, the cell index, and the
previous lane's F and h. -/
def vecRest : Int → Int → Int → Int → List Int → List (Int × Int × Bool)
  | _, _, _, _, [] => []
  | carry, b, fp, hp, x :: t =>
    let F := carry - b - 6
    let h := max x F
    (F, h, decide (fp + 6 > hp)) :: vecRest (max carry (x + b)) (b + 1) F h t

def vecRow (NEG : Int) : List Int → List (Int × Int × Bool)
  | [] => []
  | x :: t => (NEG, max x NEG, false) :: vecRest (max (NEG + 6) x) 1 NEG (max x NEG) t

/-- The invariant: the deletion value the scalar loop is about to form
equals carry - b - 6. -/
theorem rest_eq : ∀ (h0 : List Int) (f hl carry b : Int),
    max (hl - 7) (f - 1) = carry - b - 6 →
    scalarRest f hl h0 = vecRest carry b f hl h0
  | [], _, _, _, _, _ => rfl
  | x :: t, f, hl, carry, b, hinv => by
    simp only [scalarRest, vecRest]
    have hbit : decide (f - 1 > hl - 7) = decide (f + 6 > hl) := by
      by_cases h : f + 6 > hl
      · have h' : f - 1 > hl - 7 := by omega
        simp [h, h']
      · have h' : ¬ (f - 1 > hl - 7) := by omega
        simp [h, h']
    have hrest := rest_eq t (carry - b - 6) (max x (carry - b - 6)) (max carry (x + b)) (b + 1)
      (by omega)
    rw [hinv, hbit, hrest]

/-- **The vector fill equals the scalar fill** in every cell of a row of
any width, for any integer inputs. -/
theorem vec_eq_scalar (NEG : Int) (h0 : List Int) : scalarRow NEG h0 = vecRow NEG h0 := by
  cases h0 with
  | nil => rfl
  | cons x t =>
    simp only [scalarRow, vecRow]
    rw [rest_eq t NEG (max x NEG) (max (NEG + 6) x) 1 (by omega)]

end Certa.AffineFill
