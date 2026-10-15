import Certa.LemmaE

/-
Lemma W (THEORY.md §4): an alignment with gap total G ≥ 1 loses at least
6 + G, so one scoring above b has G ≤ pad(b) = max(0, L − b − 7). Every
certificate's band (diagonals ± pad) rests on it.

On Lemma E's alignments: c₁ and c₂ clipped bases around a transcript; the
gap total counts inserted and deleted bases.
-/
namespace Certa

def stGap : St → Nat
  | .m => 0
  | .x => 0
  | .i => 1

/-- The gap total: inserted and deleted bases. -/
def gaps : List Base → Nat
  | [] => 0
  | b :: t => stGap b.s + b.del + gaps t

/-- 6 for opening a gap, unless an insertion run is already open (`r`). -/
def openW (r : Bool) (G : Nat) : Nat := if r then 0 else if G = 0 then 0 else 6

theorem gap_step (r : Bool) (b : Base) (gt : Nat) :
    stGap b.s + b.del + openW r (stGap b.s + b.del + gt) ≤
      stLoss r b.s + delLoss b.del + openW (cont b) gt := by
  obtain ⟨s, d⟩ := b
  cases s <;> cases d <;> cases r <;>
    simp [openW, stGap, stLoss, delLoss, cont_m, cont_x, cont_i0, cont_i_succ] <;>
    (repeat' split) <;> omega

theorem loss_ge_gaps (r : Bool) (bs : List Base) : gaps bs + openW r (gaps bs) ≤ lossAux r bs := by
  induction bs generalizing r with
  | nil => cases r <;> simp [gaps, openW, lossAux]
  | cons b t ih =>
    have h1 := ih (cont b)
    have h2 := gap_step r b (gaps t)
    simp only [gaps, lossAux]
    omega

/-- **Lemma W.** An alignment scoring above b has gap total at most
max(0, L − b − 7), L its read length. -/
theorem lemmaW (c₁ c₂ : Nat) (bs : List Base) (b : Int) (hs : localScore c₁ c₂ bs > b) :
    (gaps bs : Int) ≤ max 0 ((c₁ + c₂ + span bs : Nat) - b - 7) := by
  have hspan := span_eq false bs
  have hl := loss_ge_gaps false bs
  simp only [openW] at hl
  simp only [localScore, clipPen] at hs
  by_cases hg : gaps bs = 0
  · rw [hg]; omega
  · simp [hg] at hl
    have : ((gaps bs : Nat) : Int) ≤ ((c₁ + c₂ + span bs : Nat) : Int) - b - 7 := by
      split at hs <;> split at hs <;> omega
    exact Int.le_trans this (Int.le_max_right _ _)

end Certa
