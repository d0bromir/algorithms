import Certa.LemmaL

/-
Lemma E (THEORY.md §4; checked on every clipped transcript of 8-base reads
by tests/verify_lemmas.cpp): a local alignment whose end-to-end completion
has d ≥ 1 edits scores at most edit_bound(L, d) = L − d − (d ≥ 2 ? 5 : 4),
and this is attained. Here for reads of any length.

The alignment: c₁ and c₂ clipped read bases, and between them a transcript
in Lemma L's model. Its score is BWA-MEM's: +1 a match, −4 a mismatch,
−(6 + k) a gap of k (a run of k inserted bases, or k deleted ones), −5 per
clipped end. Its completion aligns the clipped bases too, so it has at most
D = mismatches + inserted + deleted + clipped bases edits (every clipped base
a mismatch, in the worst case); edit_bound decreases in d.
-/
namespace Certa

def stScore (r : Bool) : St → Int
  | .m => 1
  | .x => -4
  | .i => if r then -1 else -7

def delScore (k : Nat) : Int := if k = 0 then 0 else -(6 + (k : Int))

def scoreAux : Bool → List Base → Int
  | _, [] => 0
  | r, b :: t => stScore r b.s + delScore b.del + scoreAux (cont b) t

def span : List Base → Nat
  | [] => 0
  | _ :: t => 1 + span t

def stEdits : St → Nat
  | .m => 0
  | .x => 1
  | .i => 1

/-- The edits of a transcript: its mismatches, inserted and deleted bases. -/
def edits : List Base → Nat
  | [] => 0
  | b :: t => stEdits b.s + b.del + edits t

/-- The loss is the span minus the score. -/
theorem span_eq (r : Bool) (bs : List Base) : (span bs : Int) = scoreAux r bs + lossAux r bs := by
  induction bs generalizing r with
  | nil => simp [span, scoreAux, lossAux]
  | cons b t ih =>
    have := ih (cont b)
    obtain ⟨s, d⟩ := b
    by_cases hd : d = 0 <;> cases s <;> cases r <;>
      simp [span, scoreAux, lossAux, stScore, stLoss, delScore, delLoss, hd] at this ⊢ <;> omega

/-- A clipped end costs 5. -/
def clipPen (c : Nat) : Nat := if c = 0 then 0 else 5

/-- The loss beyond the edits: at least 4 for one edit, 5 for two or more. -/
def bonus : Nat → Nat
  | 0 => 0
  | 1 => 4
  | _ + 2 => 5

/-- Inside an insertion run that started before the list, each further
inserted base adds only 1 beyond its edit: the bound there is min(e, 2). -/
def bonusT : Nat → Nat
  | 0 => 0
  | 1 => 1
  | _ + 2 => 2

def bonusR : Bool → Nat → Nat
  | true, e => bonusT e
  | false, e => bonus e

theorem bonus_spec (e : Nat) :
    (e = 0 ∧ bonus e = 0) ∨ (e = 1 ∧ bonus e = 4) ∨ (2 ≤ e ∧ bonus e = 5) := by
  match e with
  | 0 => simp [bonus]
  | 1 => simp [bonus]
  | k + 2 => simp [bonus]

theorem bonusT_spec (e : Nat) :
    (e = 0 ∧ bonusT e = 0) ∨ (e = 1 ∧ bonusT e = 1) ∨ (2 ≤ e ∧ bonusT e = 2) := by
  match e with
  | 0 => simp [bonusT]
  | 1 => simp [bonusT]
  | k + 2 => simp [bonusT]

theorem cont_m (d : Nat) : cont ⟨.m, d⟩ = false := rfl
theorem cont_x (d : Nat) : cont ⟨.x, d⟩ = false := rfl
theorem cont_i0 : cont ⟨.i, 0⟩ = true := rfl
theorem cont_i_succ (k : Nat) : cont ⟨.i, k + 1⟩ = false := rfl

theorem loss_step (r : Bool) (b : Base) (et : Nat) :
    stEdits b.s + b.del + bonusR r (stEdits b.s + b.del + et) ≤
      stLoss r b.s + delLoss b.del + bonusR (cont b) et := by
  obtain ⟨s, d⟩ := b
  have h1 := bonus_spec (stEdits s + d + et)
  have h2 := bonusT_spec (stEdits s + d + et)
  have h3 := bonus_spec et
  have h4 := bonusT_spec et
  cases s <;> cases d <;> cases r <;>
    simp [bonusR, stEdits, stLoss, delLoss, cont_m, cont_x, cont_i0, cont_i_succ] at h1 h2 h3 h4 ⊢ <;>
    omega

theorem loss_ge (r : Bool) (bs : List Base) : edits bs + bonusR r (edits bs) ≤ lossAux r bs := by
  induction bs generalizing r with
  | nil => cases r <;> simp [edits, bonusR, bonus, bonusT, lossAux]
  | cons b t ih =>
    have h1 := ih (cont b)
    have h2 := loss_step r b (edits t)
    simp only [edits, lossAux]
    omega

/-- c clipped bases lose c + 5 against c edits: at least the bonus. -/
theorem bonus_le_clip (c : Nat) : bonus c ≤ clipPen c := by
  have := bonus_spec c
  by_cases h : c = 0 <;> simp [clipPen, h] at this ⊢ <;> omega

theorem bonus_add (a b : Nat) : bonus (a + b) ≤ bonus a + bonus b := by
  have := bonus_spec a
  have := bonus_spec b
  have := bonus_spec (a + b)
  omega

/-- BWA-MEM's bound from an edit count (certa/core.h edit_bound). -/
def editBound (L d : Nat) : Int := (L : Int) - d - (if 2 ≤ d then 5 else 4)

/-- The score of a local alignment: c₁ and c₂ clipped bases around `bs`. -/
def localScore (c₁ c₂ : Nat) (bs : List Base) : Int :=
  scoreAux false bs - clipPen c₁ - clipPen c₂

theorem editBound_eq (L d : Nat) (h : 1 ≤ d) : editBound L d = (L : Int) - d - bonus d := by
  have := bonus_spec d
  by_cases h2 : 2 ≤ d <;> simp [editBound, h2] <;> omega

/-- **Lemma E**, for D = the transcript's edits plus the clipped bases. -/
theorem lemmaE_D (c₁ c₂ : Nat) (bs : List Base) (hD : 1 ≤ edits bs + c₁ + c₂) :
    localScore c₁ c₂ bs ≤ editBound (c₁ + c₂ + span bs) (edits bs + c₁ + c₂) := by
  have hs := span_eq false bs
  have hl : edits bs + bonus (edits bs) ≤ lossAux false bs := by simpa [bonusR] using loss_ge false bs
  have h1 := bonus_le_clip c₁
  have h2 := bonus_le_clip c₂
  have hb1 := bonus_add (edits bs) c₁
  have hb2 := bonus_add (edits bs + c₁) c₂
  -- the loss covers the edits and the bonus of all of them
  have core : edits bs + bonus (edits bs + c₁ + c₂) ≤ lossAux false bs + clipPen c₁ + clipPen c₂ := by
    omega
  rw [editBound_eq _ _ hD]
  simp only [localScore]
  omega

theorem editBound_anti (L d d' : Nat) (h1 : 1 ≤ d) (h : d ≤ d') : editBound L d' ≤ editBound L d := by
  simp only [editBound]
  by_cases a : 2 ≤ d <;> by_cases b : 2 ≤ d' <;> simp [a, b] <;> omega

/-- **Lemma E.** If the completion has d ≥ 1 edits, at most D (the clipped
bases counted as mismatches), the local score is at most edit_bound(L, d). -/
theorem lemmaE (c₁ c₂ : Nat) (bs : List Base) (d : Nat) (h1 : 1 ≤ d) (hd : d ≤ edits bs + c₁ + c₂) :
    localScore c₁ c₂ bs ≤ editBound (c₁ + c₂ + span bs) d :=
  Int.le_trans (lemmaE_D c₁ c₂ bs (Nat.le_trans h1 hd))
    (editBound_anti _ d _ h1 hd)

theorem matches_span (k : Nat) : span (List.replicate k (⟨.m, 0⟩ : Base)) = k := by
  induction k with
  | zero => rfl
  | succ k ih => simp only [List.replicate_succ, span, ih]; omega

theorem matches_edits (k : Nat) : edits (List.replicate k (⟨.m, 0⟩ : Base)) = 0 := by
  induction k with
  | zero => rfl
  | succ k ih => simp only [List.replicate_succ, edits, stEdits, ih]

theorem matches_score (k : Nat) : scoreAux false (List.replicate k (⟨.m, 0⟩ : Base)) = k := by
  induction k with
  | zero => rfl
  | succ k ih =>
    simp [List.replicate_succ, scoreAux, cont_m, stScore, delScore, ih]
    omega

/-- The bound is attained, by an alignment of n ≥ 1 aligned bases: for d = 1
one of them a mismatch, for d ≥ 2 all matches and d clipped bases. -/
theorem lemmaE_attained (n d : Nat) (hn : 1 ≤ n) (h1 : 1 ≤ d) :
    ∃ c₁ c₂ bs, 1 ≤ span bs ∧ edits bs + c₁ + c₂ = d ∧ c₁ + c₂ + span bs = n + d ∧
      localScore c₁ c₂ bs = editBound (n + d) d := by
  by_cases hd1 : d = 1
  · subst hd1
    refine ⟨0, 0, ⟨.x, 0⟩ :: List.replicate n ⟨.m, 0⟩, ?_, ?_, ?_, ?_⟩
    · simp only [span, matches_span]; omega
    · simp only [edits, stEdits, matches_edits]
    · simp only [span, matches_span]; omega
    · (simp [localScore, scoreAux, cont_x, stScore, delScore, matches_score, clipPen, editBound]) <;> omega
  · have h2 : 2 ≤ d := by omega
    have hd0 : d ≠ 0 := by omega
    refine ⟨d, 0, List.replicate n ⟨.m, 0⟩, ?_, ?_, ?_, ?_⟩
    · simp only [matches_span]; omega
    · simp only [matches_edits]; omega
    · simp only [matches_span]; omega
    · (simp [localScore, matches_score, clipPen, editBound, hd0, h2]) <;> omega

end Certa
