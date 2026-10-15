import Certa.LemmaL

/-
Lemma H (THEORY.md §11 list; checked on every layout of 10-base reads and on
three parts of 4 in 12-base reads by tests/verify_lemmas.cpp): with parts of
length ≥ 4 and one-mismatch probes, every part that no probe witnesses costs
at least 7, and these costs are disjoint. Here for transcripts of any length.

A part is witnessed when an ungapped occurrence with at most one mismatch
lies on the alignment's path, read off the transcript as in verify_lemmas
(witnessed1): no deletion inside the part, and its bases all matches except
at most one mismatch, or except one inserted base at an end of the part (its
other bases then lie on one diagonal, with no mismatch).

The proof reuses Lemma L's credits (`credits_le_loss`): every unwitnessed
part of length ≥ 4 collects at least 7 (`part_credit7`).
-/
namespace Certa

/-- No deletion between two bases of the part. -/
def noInnerDel : List Base → Bool
  | [] => true
  | [_] => true
  | b :: c :: t => b.del == 0 && noInnerDel (c :: t)

def nX : List Base → Nat
  | [] => 0
  | b :: t => (if b.s = .x then 1 else 0) + nX t

def nI : List Base → Nat
  | [] => 0
  | b :: t => (if b.s = .i then 1 else 0) + nI t

def headSt : List Base → Option St
  | [] => none
  | b :: _ => some b.s

def lastSt : List Base → Option St
  | [] => none
  | [b] => some b.s
  | _ :: c :: t => lastSt (c :: t)

/-- Witnessed by a one-mismatch probe (verify_lemmas.cpp `witnessed1`). -/
def Witnessed (p : List Base) : Prop :=
  noInnerDel p = true ∧
    (nX p + nI p = 0 ∨ (nI p = 0 ∧ nX p = 1) ∨
      (nI p = 1 ∧ nX p = 0 ∧ (headSt p = some .i ∨ lastSt p = some .i)))

instance (p : List Base) : Decidable (Witnessed p) := by unfold Witnessed; infer_instance

/-- An exact part has no deletion inside, no mismatch and no inserted base. -/
theorem exact_facts (l : List Base) (h : exact l = true) :
    noInnerDel l = true ∧ nX l = 0 ∧ nI l = 0 := by
  induction l with
  | nil => simp [noInnerDel, nX, nI]
  | cons b t ih =>
    cases t with
    | nil =>
      have hb : b.s = .m := by simpa [exact] using h
      simp [noInnerDel, nX, nI, hb]
    | cons c t' =>
      simp [exact] at h
      obtain ⟨⟨hs, hd⟩, ht⟩ := h
      have := ih ht
      simp [noInnerDel, nX, nI, hs, hd] at this ⊢
      exact this

theorem lastSt_cons (b : Base) (l : List Base) (h : l ≠ []) : lastSt (b :: l) = lastSt l := by
  cases l with
  | nil => exact absurd rfl h
  | cons c t => rfl

theorem credit_open_close (b : Base) (nx : Option St) (h : b.s = .i) (hc : ¬(b.del = 0 ∧ nx = some .i)) :
    8 ≤ credit false b nx := by
  simp only [credit, h, hc, Bool.false_eq_true, ↓reduceIte]
  omega

/-- A part whose bases after a match `b` are witnessed, while the whole part
is not, starts with a lone inserted base followed by a non-inserted one. -/
theorem witnessed_after_match (b c : Base) (t : List Base) (hb : b.s = .m) (hd : b.del = 0)
    (hw : Witnessed (c :: t)) (hnw : ¬Witnessed (b :: c :: t)) :
    c.s = .i ∧ c.del = 0 ∧ ∃ e t'', t = e :: t'' ∧ e.s ≠ .i := by
  have hnd : noInnerDel (b :: c :: t) = noInnerDel (c :: t) := by simp [noInnerDel, hd]
  have hx : nX (b :: c :: t) = nX (c :: t) := by simp [nX, hb]
  have hi : nI (b :: c :: t) = nI (c :: t) := by simp [nI, hb]
  have hl : lastSt (b :: c :: t) = lastSt (c :: t) := lastSt_cons b (c :: t) (by simp)
  have hh : headSt (b :: c :: t) = some .m := by simp [headSt, hb]
  unfold Witnessed at hw hnw
  rw [hnd, hx, hi, hl, hh] at hnw
  obtain ⟨hw1, hw2⟩ := hw
  -- the only way out: one inserted base, at the head of c :: t but not at its end
  have hci : c.s = .i ∧ lastSt (c :: t) ≠ some .i ∧ nI (c :: t) = 1 := by
    rcases hw2 with h0 | ⟨h1, h2⟩ | ⟨h1, h2, h3 | h3⟩
    · exact absurd ⟨hw1, Or.inl h0⟩ hnw
    · exact absurd ⟨hw1, Or.inr (Or.inl ⟨h1, h2⟩)⟩ hnw
    · refine ⟨by simpa [headSt] using h3, fun hl' => hnw ⟨hw1, Or.inr (Or.inr ⟨h1, h2, Or.inr hl'⟩)⟩, h1⟩
    · exact absurd ⟨hw1, Or.inr (Or.inr ⟨h1, h2, Or.inr h3⟩)⟩ hnw
  obtain ⟨hcs, hlast, hn1⟩ := hci
  cases t with
  | nil => simp [lastSt, hcs] at hlast
  | cons e t'' =>
    refine ⟨hcs, ?_, e, t'', rfl, ?_⟩
    · simp [noInnerDel] at hw1; exact hw1.1
    · intro he
      simp [nI, hcs, he] at hn1

/-- Every unwitnessed part collects at least 7 credits, if no run enters it
from the left (`r = false`) or it has at least 4 bases. -/
theorem part_credit7 (r : Bool) (seg : List Base) (nx : Option St) (hw : ¬Witnessed seg)
    (hlen : r = false ∨ 4 ≤ seg.length) : 7 ≤ sum (credits r seg nx) := by
  induction seg generalizing r with
  | nil => exact absurd (by decide) hw
  | cons b t ih =>
    rw [sum_cons]
    cases t with
    | nil =>
      exfalso; apply hw
      cases hs : b.s <;> simp [Witnessed, noInnerDel, nX, nI, headSt, hs]
    | cons c t' =>
      have hn : nextSt (c :: t') nx = some c.s := rfl
      rw [hn]
      by_cases hd : b.del = 0
      case neg =>
        have := credit_del r b (some c.s)
        have := delLoss_pos b.del hd
        omega
      cases hs : b.s with
      | x =>
        have h5 := credit_x r b (some c.s) hs
        by_cases he : exact (c :: t') = true
        · exfalso; apply hw
          obtain ⟨h1, h2, h3⟩ := exact_facts _ he
          refine ⟨by simp [noInnerDel, hd, h1], Or.inr (Or.inl ⟨?_, ?_⟩)⟩
          · simp [nI, hs] at h3 ⊢; exact h3
          · simp [nX, hs] at h2 ⊢; exact h2
        · have hcb : cont b = false := by simp [cont, hs]
          have := part_credit false (c :: t') nx (by simpa using he) (Or.inl rfl)
          rw [hcb]
          omega
      | m =>
        have hcb : cont b = false := by simp [cont, hs]
        rw [hcb]
        by_cases hwt : Witnessed (c :: t')
        · obtain ⟨hci, hcd, e, t'', rfl, he⟩ := witnessed_after_match b c t' hs hd hwt hw
          rw [sum_cons]
          have hn2 : nextSt (e :: t'') nx = some e.s := rfl
          rw [hn2]
          have := credit_open_close c (some e.s) hci (by simp [he])
          omega
        · have := ih false hwt (Or.inl rfl)
          omega
      | i =>
        cases r with
        | false =>
          by_cases hci : c.s = .i
          · have h5 := credit_open b (some c.s) hs
            rw [sum_cons]
            have := credit_i (cont b) c (nextSt t' nx) hci
            omega
          · have := credit_open_close b (some c.s) hs (by simp [hci])
            omega
        | true =>
          have hlen4 : 4 ≤ (b :: c :: t').length := by simpa using hlen
          have hcb : cont b = true := by simp [cont, hs, hd]
          rw [hcb]
          by_cases hci : c.s = .i
          · -- the run goes on into c
            have h2 := credit_i true b (some c.s) hs
            cases t' with
            | nil => simp at hlen4
            | cons d t'' =>
              rw [sum_cons]
              have hn2 : nextSt (d :: t'') nx = some d.s := rfl
              rw [hn2]
              by_cases hcl : c.del = 0 ∧ some d.s = some St.i
              · have hc2 := credit_i true c (some d.s) hci
                have hds : d.s = .i := by simpa using hcl.2
                have hcc : cont c = true := by simp [cont, hci, hcl.1]
                rw [hcc]
                cases t'' with
                | nil => simp at hlen4
                | cons e t3 =>
                  rw [sum_cons]
                  have hn3 : nextSt (e :: t3) nx = some e.s := rfl
                  rw [hn3]
                  by_cases hdl : d.del = 0 ∧ some e.s = some St.i
                  · have hd2 := credit_i true d (some e.s) hds
                    have hes : e.s = .i := by simpa using hdl.2
                    rw [sum_cons]
                    have := credit_i (cont d) e (nextSt t3 nx) hes
                    omega
                  · have := credit_close true d (some e.s) hds hdl
                    omega
              · have := credit_close true c (some d.s) hci hcl
                omega
          · -- b closes its run; the rest is inexact, or the part is witnessed
            have h5 := credit_close true b (some c.s) hs (by simp [hci])
            by_cases he : exact (c :: t') = true
            · exfalso; apply hw
              obtain ⟨h1, h2, h3⟩ := exact_facts _ he
              refine ⟨by simp [noInnerDel, hd, h1], Or.inr (Or.inr ⟨?_, ?_, Or.inl (by simp [headSt, hs])⟩)⟩
              · simp [nI, hs] at h3 ⊢; exact h3
              · simp [nX, hs] at h2 ⊢; exact h2
            · have := part_credit true (c :: t') nx (by simpa using he) (Or.inr (by simp at hlen4 ⊢; omega))
              omega

theorem parts_credit7 (segs : List (List Base × List Base)) (tail : List Base) (r : Bool) (nx : Option St)
    (h : ∀ gp ∈ segs, 4 ≤ gp.2.length) :
    7 * (segs.filter (fun gp => !decide (Witnessed gp.2))).length ≤ sum (credits r (assemble segs tail) nx) := by
  induction segs generalizing r with
  | nil => simp
  | cons gp rest ih =>
    obtain ⟨g, p⟩ := gp
    have hp : 4 ≤ p.length := h (g, p) (by simp)
    have hrest : ∀ gp ∈ rest, 4 ≤ gp.2.length := fun gp hm => h gp (by simp [hm])
    simp only [assemble, List.append_assoc, credits_append, sum_append]
    have ihr := ih (flagAfter (flagAfter r g) p) hrest
    by_cases hw : Witnessed p
    · simp [List.filter, hw]
      omega
    · have := part_credit7 (flagAfter r g) p (nextSt (assemble rest tail) nx) hw (Or.inr hp)
      simp [List.filter, hw]
      omega

/-- **Lemma H.** In a read cut into gaps and parts of length ≥ 4, the loss of
an alignment is at least 7 per part that no one-mismatch probe witnesses. -/
theorem lemmaH (segs : List (List Base × List Base)) (tail : List Base)
    (h : ∀ gp ∈ segs, 4 ≤ gp.2.length) :
    7 * (segs.filter (fun gp => !decide (Witnessed gp.2))).length ≤ loss (assemble segs tail) :=
  Nat.le_trans (parts_credit7 segs tail false none h) (credits_le_loss _ none)

/-- Parts of length 3 do not suffice: seven inserted bases over three parts of
3 (a match, then the run, then a match) have loss 20 < 21. -/
theorem lemmaH_length_three :
    ∃ segs tail, (∀ gp ∈ segs, 3 ≤ gp.2.length) ∧
      loss (assemble segs tail) < 7 * (segs.filter (fun gp => !decide (Witnessed gp.2))).length := by
  refine ⟨[([], [⟨.m, 0⟩, ⟨.i, 0⟩, ⟨.i, 0⟩]), ([], [⟨.i, 0⟩, ⟨.i, 0⟩, ⟨.i, 0⟩]),
    ([], [⟨.i, 0⟩, ⟨.i, 0⟩, ⟨.m, 0⟩])], [], ?_, ?_⟩
  · decide
  · decide

end Certa
