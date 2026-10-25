/-
Lemma L (THEORY.md §4; checked on every transcript of up to 9 read bases by
tests/verify_lemmas.cpp): in an alignment, every part of the read of length
≥ 3 that does not occur exactly costs at least 5 points, and these costs are
disjoint. Here for transcripts of any length.

The model, as in verify_lemmas.cpp. Each read base of the aligned span is a
match, a mismatch or inserted, and a deletion of `del` reference bases may
follow it (0: none). The loss (span minus score) is: a mismatch 5; an
inserted base 2, plus 6 for each run of inserted bases (consecutive, with no
deletion between them); a deletion of k ≥ 1 bases 6 + k. A part occurs
exactly when all its bases are matches and no deletion lies between two of
them.

The proof gives every base a credit: its own loss, with each run's opening
cost of 6 split as 3 to the run's first base and 3 to its last. The credits
sum to at most the loss (`credits_le_loss`); every inexact part of length
≥ 3 collects at least 5 (`part_credit`); parts are disjoint (`lemmaL`).
`lemmaL_length_two` shows that parts of length 2 do not suffice.
-/
namespace Certa

inductive St where
  | m  -- match
  | x  -- mismatch
  | i  -- inserted
  deriving DecidableEq, Repr

structure Base where
  s : St
  del : Nat
  deriving DecidableEq, Repr

def delLoss (k : Nat) : Nat := if k = 0 then 0 else 6 + k

/-- Whether an inserted base after `b` continues `b`'s run. -/
def cont (b : Base) : Bool := b.s == .i && b.del == 0

/-- `r`: the previous base continues into this one (`cont`). -/
def stLoss (r : Bool) : St → Nat
  | .m => 0
  | .x => 5
  | .i => if r then 2 else 8

def lossAux : Bool → List Base → Nat
  | _, [] => 0
  | r, b :: t => stLoss r b.s + delLoss b.del + lossAux (cont b) t

/-- The loss of a transcript: its span minus its score. -/
def loss (bs : List Base) : Nat := lossAux false bs

/-- The state of the base after `t`, given the one after the whole list. -/
def nextSt (t : List Base) (nx : Option St) : Option St :=
  match t with
  | [] => nx
  | c :: _ => some c.s

/-- A base's credit: its own loss, and for an inserted base 3 when it opens
a run and 3 when it closes one. `nx`: the next base's state. -/
def credit (r : Bool) (b : Base) (nx : Option St) : Nat :=
  delLoss b.del +
    match b.s with
    | .m => 0
    | .x => 5
    | .i => 2 + (if r then 0 else 3) + (if b.del = 0 ∧ nx = some .i then 0 else 3)

def credits : Bool → List Base → Option St → List Nat
  | _, [], _ => []
  | r, b :: t, nx => credit r b (nextSt t nx) :: credits (cont b) t nx

def sum : List Nat → Nat
  | [] => 0
  | a :: t => a + sum t

theorem sum_append (u v : List Nat) : sum (u ++ v) = sum u + sum v := by
  induction u with
  | nil => simp [sum]
  | cons a t ih => simp [sum, ih]; omega

theorem sum_cons (r : Bool) (b : Base) (t : List Base) (nx : Option St) :
    sum (credits r (b :: t) nx) = credit r b (nextSt t nx) + sum (credits (cont b) t nx) := rfl

/-- One base: its credit, plus 3 owed to the next base when the run goes on
into it, is at most its loss, plus 3 when a run enters it (`r`). -/
theorem step (r : Bool) (b : Base) (nxt : Option St) :
    credit r b nxt + (if cont b = true ∧ nxt = some .i then 3 else 0) ≤
      stLoss r b.s + delLoss b.del + (if r = true ∧ b.s = .i then 3 else 0) := by
  obtain ⟨s, d⟩ := b
  by_cases hd : d = 0 <;> by_cases hn : nxt = some .i <;> cases s <;> cases r <;>
    simp [credit, stLoss, cont, delLoss, hd, hn] <;> omega

/-- The credits never exceed the loss; a run open before the list (r) may
give 3 more, to its closing base. -/
theorem credits_le (r : Bool) (bs : List Base) (nx : Option St) :
    sum (credits r bs nx) ≤
      lossAux r bs + (if r = true ∧ (bs.head?.map (·.s)) = some .i then 3 else 0) := by
  induction bs generalizing r with
  | nil => simp [credits, sum, lossAux]
  | cons b t ih =>
    have hs := step r b (nextSt t nx)
    have hb : (if r = true ∧ ((b :: t).head?.map (·.s)) = some .i then 3 else 0) =
        (if r = true ∧ b.s = .i then 3 else 0) := by simp
    rw [sum_cons, hb]
    simp only [lossAux]
    cases t with
    | nil =>
      have : sum (credits (cont b) [] nx) = 0 := rfl
      have : lossAux (cont b) [] = 0 := rfl
      omega
    | cons c t' =>
      have h := ih (cont b)
      have hc : ((c :: t').head?.map (·.s)) = some c.s := rfl
      have hn : nextSt (c :: t') nx = some c.s := rfl
      rw [hc] at h
      rw [hn] at hs ⊢
      omega

theorem credits_le_loss (bs : List Base) (nx : Option St) : sum (credits false bs nx) ≤ loss bs := by
  have := credits_le false bs nx
  simpa [loss] using this

/-- The flag after a list: whether its last base continues into the next. -/
def flagAfter : Bool → List Base → Bool
  | r, [] => r
  | _, b :: t => flagAfter (cont b) t

theorem nextSt_append (t v : List Base) (nx : Option St) : nextSt (t ++ v) nx = nextSt t (nextSt v nx) := by
  cases t <;> rfl

theorem credits_append (r : Bool) (u v : List Base) (nx : Option St) :
    credits r (u ++ v) nx = credits r u (nextSt v nx) ++ credits (flagAfter r u) v nx := by
  induction u generalizing r with
  | nil => rfl
  | cons b t ih => simp [credits, flagAfter, ih, nextSt_append]

/-- A part occurs exactly: all its bases are matches, and no deletion lies
between two of them (the last base's deletion is outside the part). -/
def exact : List Base → Bool
  | [] => true
  | [b] => b.s == .m
  | b :: c :: t => b.s == .m && b.del == 0 && exact (c :: t)

theorem delLoss_pos (d : Nat) (h : d ≠ 0) : 7 ≤ delLoss d := by simp [delLoss, h]; omega

theorem credit_del (r : Bool) (b : Base) (nx : Option St) : delLoss b.del ≤ credit r b nx :=
  Nat.le_add_right _ _

theorem credit_x (r : Bool) (b : Base) (nx : Option St) (h : b.s = .x) : 5 ≤ credit r b nx := by
  simp [credit, h]

theorem credit_i (r : Bool) (b : Base) (nx : Option St) (h : b.s = .i) : 2 ≤ credit r b nx := by
  simp [credit, h]; omega

theorem credit_open (b : Base) (nx : Option St) (h : b.s = .i) : 5 ≤ credit false b nx := by
  simp [credit, h]; omega

theorem credit_close (r : Bool) (b : Base) (nx : Option St) (h : b.s = .i)
    (hc : ¬(b.del = 0 ∧ nx = some .i)) : 5 ≤ credit r b nx := by
  simp only [credit, h, hc, ↓reduceIte]
  omega

/-- Every inexact part collects at least 5 credits, if no run enters it from
the left (`r = false`) or it has at least 3 bases. -/
theorem part_credit (r : Bool) (seg : List Base) (nx : Option St) (hx : exact seg = false)
    (hlen : r = false ∨ 3 ≤ seg.length) : 5 ≤ sum (credits r seg nx) := by
  induction seg generalizing r with
  | nil => simp [exact] at hx
  | cons b t ih =>
    rw [sum_cons]
    cases t with
    | nil =>
      have hr : r = false := by simpa using hlen
      subst hr
      cases hs : b.s with
      | m => simp [exact, hs] at hx
      | x => have := credit_x false b (nextSt [] nx) hs; omega
      | i => have := credit_open b (nextSt [] nx) hs; omega
    | cons c t' =>
      have hn : nextSt (c :: t') nx = some c.s := rfl
      rw [hn]
      cases hs : b.s with
      | x => have := credit_x r b (some c.s) hs; omega
      | m =>
        by_cases hd : b.del = 0
        · have hx' : exact (c :: t') = false := by simpa [exact, hs, hd] using hx
          have hcb : cont b = false := by simp [cont, hs]
          have := ih false hx' (Or.inl rfl)
          rw [hcb]
          omega
        · have := credit_del r b (some c.s)
          have := delLoss_pos b.del hd
          omega
      | i =>
        by_cases hd : b.del = 0
        · cases r with
          | false => have := credit_open b (some c.s) hs; omega
          | true =>
            by_cases hci : c.s = .i
            · -- the run goes on into c; a third base exists
              have hlen3 : 3 ≤ (b :: c :: t').length := by simpa using hlen
              cases t' with
              | nil => simp at hlen3
              | cons d t'' =>
                have hb2 := credit_i true b (some c.s) hs
                have hcb : cont b = true := by simp [cont, hs, hd]
                rw [hcb, sum_cons]
                have hn2 : nextSt (d :: t'') nx = some d.s := rfl
                rw [hn2]
                by_cases hcl : c.del = 0 ∧ some d.s = some St.i
                · -- the run goes on into d too: 2 + 2 + 2
                  have hc2 := credit_i true c (some d.s) hci
                  have hds : d.s = .i := by simpa using hcl.2
                  rw [sum_cons]
                  have hd2 := credit_i (cont c) d (nextSt t'' nx) hds
                  omega
                · have hc5 := credit_close true c (some d.s) hci hcl
                  omega
            · have := credit_close true b (some c.s) hs (by simp [hci])
              omega
        · have := credit_del r b (some c.s)
          have := delLoss_pos b.del hd
          omega

/-- A read cut into gaps and parts: `(gap, part)` pairs, then a tail. -/
def assemble : List (List Base × List Base) → List Base → List Base
  | [], tail => tail
  | (g, p) :: rest, tail => g ++ p ++ assemble rest tail

theorem parts_credit (segs : List (List Base × List Base)) (tail : List Base) (r : Bool) (nx : Option St)
    (h : ∀ gp ∈ segs, 3 ≤ gp.2.length) :
    5 * (segs.filter (fun gp => !exact gp.2)).length ≤ sum (credits r (assemble segs tail) nx) := by
  induction segs generalizing r with
  | nil => simp
  | cons gp rest ih =>
    obtain ⟨g, p⟩ := gp
    have hp : 3 ≤ p.length := h (g, p) (by simp)
    have hrest : ∀ gp ∈ rest, 3 ≤ gp.2.length := fun gp hm => h gp (by simp [hm])
    simp only [assemble, List.append_assoc, credits_append, sum_append]
    have ihr := ih (flagAfter (flagAfter r g) p) hrest
    by_cases hx : exact p = false
    · have := part_credit (flagAfter r g) p (nextSt (assemble rest tail) nx) hx (Or.inr hp)
      simp [List.filter, hx]
      omega
    · have hx' : exact p = true := by simpa using hx
      simp [List.filter, hx']
      omega

/-- **Lemma L.** In a read cut into gaps and parts of length ≥ 3, the loss of
an alignment is at least 5 per part that does not occur exactly. -/
theorem lemmaL (segs : List (List Base × List Base)) (tail : List Base)
    (h : ∀ gp ∈ segs, 3 ≤ gp.2.length) :
    5 * (segs.filter (fun gp => !exact gp.2)).length ≤ loss (assemble segs tail) :=
  Nat.le_trans (parts_credit segs tail false none h) (credits_le_loss _ none)

/-- Parts of length 2 do not suffice (verify_lemmas.cpp's counterexample):
a match, four inserted bases and two matches, cut into three parts of 2, have
loss 14 < 15. -/
theorem lemmaL_length_two :
    ∃ segs tail, (∀ gp ∈ segs, 2 ≤ gp.2.length) ∧
      loss (assemble segs tail) < 5 * (segs.filter (fun gp => !exact gp.2)).length := by
  refine ⟨[([], [⟨.m, 0⟩, ⟨.i, 0⟩]), ([], [⟨.i, 0⟩, ⟨.i, 0⟩]), ([], [⟨.i, 0⟩, ⟨.m, 0⟩])], [⟨.m, 0⟩, ⟨.m, 0⟩], ?_, ?_⟩
  · decide
  · decide

end Certa
