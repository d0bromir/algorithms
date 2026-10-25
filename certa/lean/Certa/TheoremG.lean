/-
Theorem G (THEORY.md §6; checked against a DP optimum on 20,000 random
instances by tests/verify_lemmas.cpp): let F be the read intervals [i, j)
occurring at most τ times in the reference. Scanning left to right and
taking at each start s the shortest interval of F starting there gives the
maximum number of disjoint intervals of F.

Here for any family F of intervals of [0, L) that is closed under extension
(`ok i j` and [i, j) ⊆ [i', j') ⊆ [0, L) give `ok i' j'`), which F_τ is: the
occurrences of a superstring are a subset of those of the string.
`greedy` is verify_lemmas.cpp's greedy (and partition_strand's, Lemma R):
from s, the least e > s with [s, e) in F, then on from e.
-/
namespace Certa.Greedy

variable (ok : Nat → Nat → Bool)

/-- The least end j, j + 1, …, j + fuel − 1 of an interval of F from s. -/
def scan (s : Nat) : Nat → Nat → Option Nat
  | 0, _ => none
  | fuel + 1, j => if ok s j then some j else scan s fuel (j + 1)

theorem scan_bounds (s : Nat) : ∀ fuel j e, scan ok s fuel j = some e → j ≤ e ∧ e < j + fuel ∧ ok s e = true
  | 0, _, _, h => by simp [scan] at h
  | fuel + 1, j, e, h => by
    simp only [scan] at h
    split at h
    · next hok => simp at h; subst h; exact ⟨Nat.le_refl _, by omega, hok⟩
    · have := scan_bounds s fuel (j + 1) e h
      exact ⟨by omega, by omega, this.2.2⟩

/-- An interval of F from s ending at k in the scanned range: the scan finds
one ending at or before k. -/
theorem scan_le (s k : Nat) (hk : ok s k = true) :
    ∀ fuel j, j ≤ k → k < j + fuel → ∃ e, scan ok s fuel j = some e ∧ e ≤ k
  | 0, _, _, h2 => by omega
  | fuel + 1, j, h1, h2 => by
    simp only [scan]
    split
    · exact ⟨j, rfl, h1⟩
    · next hj =>
      have hjk : j ≠ k := fun he => by subst he; exact hj hk
      exact scan_le s k hk fuel (j + 1) (by omega) (by omega)

/-- The greedy partition from s: the shortest interval of F at s, then on
from its end. -/
def greedy (L s : Nat) : List (Nat × Nat) :=
  match _h : scan ok s (L - s) (s + 1) with
  | none => []
  | some e => (s, e) :: greedy L e
termination_by L - s
decreasing_by
  have := scan_bounds ok s (L - s) (s + 1) e _h
  omega

theorem greedy_none (L s : Nat) (h : scan ok s (L - s) (s + 1) = none) : greedy ok L s = [] := by
  rw [greedy]; split <;> simp_all

theorem greedy_some (L s e : Nat) (h : scan ok s (L - s) (s + 1) = some e) :
    greedy ok L s = (s, e) :: greedy ok L e := by
  rw [greedy]; split <;> simp_all

/-- Disjoint intervals of F in order, the first starting at or after s, all
in [0, L). -/
def Family (L : Nat) : Nat → List (Nat × Nat) → Prop
  | _, [] => True
  | s, (i, j) :: t => s ≤ i ∧ i < j ∧ j ≤ L ∧ ok i j = true ∧ Family L j t

theorem family_mono (L : Nat) : ∀ (fam : List (Nat × Nat)) (s s' : Nat), s' ≤ s → Family ok L s fam → Family ok L s' fam
  | [], _, _, _, _ => trivial
  | (i, j) :: t, s, s', hs, ⟨h1, h2, h3, h4, h5⟩ => ⟨by omega, h2, h3, h4, h5⟩

/-- The greedy picks are disjoint intervals of F. -/
theorem greedy_family (L s : Nat) : Family ok L s (greedy ok L s) := by
  induction s using greedy.induct ok L with
  | case1 s h => rw [greedy_none ok L s h]; trivial
  | case2 s e h ih =>
    rw [greedy_some ok L s e h]
    obtain ⟨b1, b2, b3⟩ := scan_bounds ok s (L - s) (s + 1) e h
    exact ⟨Nat.le_refl _, by omega, by omega, b3, ih⟩

/-- **Theorem G.** If F is closed under extension, no disjoint family of
intervals of F starting at or after s is larger than the greedy one. -/
theorem greedy_optimal (L : Nat)
    (mono : ∀ i j i' j', ok i j = true → i' ≤ i → j ≤ j' → j' ≤ L → ok i' j' = true) :
    ∀ (fam : List (Nat × Nat)) (s : Nat), Family ok L s fam → fam.length ≤ (greedy ok L s).length
  | [], _, _ => Nat.zero_le _
  | (i, j) :: t, s, ⟨h1, h2, h3, h4, h5⟩ => by
    -- [i, j) in F extends to [s, j): the greedy interval at s ends at or before j
    have hs : ok s j = true := mono i j s j h4 h1 (Nat.le_refl _) h3
    obtain ⟨e, he, hej⟩ := scan_le ok s j hs (L - s) (s + 1) (by omega) (by omega)
    rw [greedy_some ok L s e he]
    simp only [List.length_cons]
    have := greedy_optimal L mono t e (family_mono ok L t j e hej h5)
    omega

/-- The greedy partition is a maximum family. -/
theorem greedy_max (L : Nat)
    (mono : ∀ i j i' j', ok i j = true → i' ≤ i → j ≤ j' → j' ≤ L → ok i' j' = true) :
    Family ok L 0 (greedy ok L 0) ∧
      ∀ fam, Family ok L 0 fam → fam.length ≤ (greedy ok L 0).length :=
  ⟨greedy_family ok L 0, fun fam h => greedy_optimal ok L mono fam 0 h⟩

/-! ### F_τ is closed under extension

Occurrences of a string w in a text T are the positions where w is a prefix
of what follows. A superstring u ++ w ++ v of w occurs at most as often as
w: extending to the right only removes matches at each position, and one
character to the left maps each match to the next position. -/

section Occurrences

variable {α : Type} [DecidableEq α]

def isPre : List α → List α → Bool
  | [], _ => true
  | _ :: _, [] => false
  | a :: w, b :: t => decide (a = b) && isPre w t

/-- The occurrences of w in T (w = [] occurs at all |T| + 1 positions). -/
def occ : List α → List α → Nat
  | [], w => if isPre w [] then 1 else 0
  | x :: T, w => (if isPre w (x :: T) then 1 else 0) + occ T w

theorem isPre_append (w v : List α) : ∀ X, isPre (w ++ v) X = true → isPre w X = true := by
  induction w with
  | nil => intro X _; cases X <;> rfl
  | cons a w ih =>
    intro X h
    cases X with
    | nil => simp [isPre] at h
    | cons b t =>
      simp [isPre] at h ⊢
      exact ⟨h.1, ih t h.2⟩

theorem occ_append_right (w v : List α) : ∀ T, occ T (w ++ v) ≤ occ T w := by
  intro T
  induction T with
  | nil =>
    simp only [occ]
    have := isPre_append w v []
    split <;> split <;> simp_all
  | cons x T ih =>
    simp only [occ]
    have := isPre_append w v (x :: T)
    split <;> split <;> simp_all <;> omega

theorem occ_tail (w : List α) : ∀ T, occ T.tail w ≤ occ T w
  | [] => Nat.le_refl _
  | x :: T => by simp only [List.tail, occ]; omega

theorem isPre_cons_le (c x : α) (w T : List α) :
    (if isPre (c :: w) (x :: T) then 1 else 0) ≤ (if isPre w T then 1 else 0) := by
  by_cases hc : c = x <;> simp [isPre, hc]

theorem occ_cons (c : α) (w : List α) : ∀ T, occ T (c :: w) ≤ occ T.tail w
  | [] => by simp [occ, isPre]
  | x :: T => by
    have ih := occ_cons c w T
    have hp := isPre_cons_le c x w T
    cases T with
    | nil => simp [occ, isPre]; split <;> simp_all
    | cons y T' =>
      show (if isPre (c :: w) (x :: y :: T') then 1 else 0) + occ (y :: T') (c :: w) ≤
        (if isPre w (y :: T') then 1 else 0) + occ T' w
      have ih' : occ (y :: T') (c :: w) ≤ occ T' w := ih
      omega

theorem occ_prepend (u w : List α) : ∀ T, occ T (u ++ w) ≤ occ T w := by
  induction u with
  | nil => intro T; exact Nat.le_refl _
  | cons c u ih =>
    intro T
    have h1 := occ_cons c (u ++ w) T
    have h2 := occ_tail (u ++ w) T
    have h3 := ih T
    simp only [List.cons_append]
    omega

/-- A superstring occurs at most as often as the string. -/
theorem occ_super (u w v T : List α) : occ T (u ++ w ++ v) ≤ occ T w := by
  have h1 := occ_append_right (u ++ w) v T
  have h2 := occ_prepend u w T
  omega

/-- The read interval [i, j). -/
def slice (r : List α) (i j : Nat) : List α := (r.drop i).take (j - i)

omit [DecidableEq α] in
theorem slice_split (r : List α) (i' i j j' : Nat) (h1 : i' ≤ i) (h2 : i ≤ j) (h3 : j ≤ j') :
    slice r i' j' = slice r i' i ++ slice r i j ++ slice r j j' := by
  simp only [slice]
  have e1 : j' - i' = (i - i') + (j' - i) := by omega
  have e2 : j' - i = (j - i) + (j' - j) := by omega
  rw [e1, List.take_add, List.drop_drop, e2, List.take_add, List.drop_drop, List.append_assoc]
  have e3 : i' + (i - i') = i := by omega
  have e4 : i + (j - i) = j := by omega
  simp only [e3, e4]

/-- **Theorem G for F_τ.** With T the reference and r the read, the
intervals of r occurring at most τ times in T are closed under extension,
so the greedy partition is a maximum family of them. -/
theorem greedy_max_tau (T r : List α) (τ : Nat) :
    let ok := fun i j => decide (occ T (slice r i j) ≤ τ)
    Family ok r.length 0 (greedy ok r.length 0) ∧
      ∀ fam, Family ok r.length 0 fam → fam.length ≤ (greedy ok r.length 0).length := by
  intro ok
  apply greedy_max
  intro i j i' j' h hi hj _
  simp only [ok, decide_eq_true_eq] at h ⊢
  by_cases hij : i ≤ j
  · rw [slice_split r i' i j j' hi hij hj]
    exact Nat.le_trans (occ_super _ _ _ T) h
  · -- an empty interval [i, j) with j < i: its slice is [], which occurs |T| + 1 times
    have he : slice r i j = [] := by simp [slice]; omega
    rw [he] at h
    have : occ T (slice r i' j') ≤ occ T ([] : List α) := by
      have := occ_super (slice r i' j') [] [] T
      simpa using this
    omega

end Occurrences

end Certa.Greedy
