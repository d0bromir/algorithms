/-
Lemma S, second claim (THEORY.md §4; checked on every layout of 24-base reads,
every counted set and extra ≤ 2 by tests/verify_lemmas.cpp): for parts of
length ≥ 5, Lemma L's bound over a span is largest on the full span [0, L).
It makes mutant L6 (shifted span start candidates) equivalent. Here for reads
of any length.

The bound for a span [x, y): (y − x) − 5 max(0, k − extra) − 5 [x > 0] −
5 [y < L], k the counted parts inside the span (verify_lemmas.cpp
span_value). Shrinking the span to [x, y) drops the parts that start before
x or end after y; parts of length ≥ 5 that start before x number at most
(x + 4) / 5, so they save at most x + 4 < x + 5, the bases and the clip
lost; symmetrically at the end.
-/
namespace Certa.SpanS

structure Part where
  off : Nat
  len : Nat
  counted : Bool
  deriving DecidableEq, Repr

/-- Disjoint parts in order, each at least `m` long, starting at or after
`lo` and ending at or before `L`. -/
def fits (m : Nat) : Nat → Nat → List Part → Bool
  | lo, L, [] => decide (lo ≤ L)
  | lo, L, p :: rest => decide (lo ≤ p.off) && decide (m ≤ p.len) && fits m (p.off + p.len) L rest

/-- The counted parts inside [x, y). -/
def kIn (x y : Nat) : List Part → Nat
  | [] => 0
  | p :: t => (if p.counted = true ∧ x ≤ p.off ∧ p.off + p.len ≤ y then 1 else 0) + kIn x y t

def kAll : List Part → Nat
  | [] => 0
  | p :: t => (if p.counted = true then 1 else 0) + kAll t

/-- The parts starting before x, and those ending after y. -/
def nLeft (x : Nat) : List Part → Nat
  | [] => 0
  | p :: t => (if p.off < x then 1 else 0) + nLeft x t

def nRight (y : Nat) : List Part → Nat
  | [] => 0
  | p :: t => (if y < p.off + p.len then 1 else 0) + nRight y t

/-- Lemma L's bound for the span [x, y) (verify_lemmas.cpp `span_value`). -/
def spanValue (L : Nat) (ps : List Part) (e x y : Nat) : Int :=
  (y : Int) - x - 5 * ((kIn x y ps - e : Nat) : Int) - (if 0 < x then 5 else 0) - (if y < L then 5 else 0)

theorem fits_len (m lo L : Nat) (ps : List Part) (h : fits m lo L ps = true) :
    lo + m * ps.length ≤ L := by
  induction ps generalizing lo with
  | nil => simpa [fits] using h
  | cons p t ih =>
    simp [fits] at h
    obtain ⟨⟨h1, h2⟩, h3⟩ := h
    have := ih _ h3
    simp [List.length]
    rw [Nat.mul_succ]
    have : m ≤ p.len := h2
    omega

/-- A counted part not inside [x, y) starts before x or ends after y. -/
theorem counted_split (x y : Nat) (ps : List Part) : kAll ps ≤ kIn x y ps + nLeft x ps + nRight y ps := by
  induction ps with
  | nil => simp [kAll, kIn, nLeft, nRight]
  | cons p t ih =>
    simp only [kAll, kIn, nLeft, nRight]
    cases p.counted <;> simp <;> (repeat' split) <;> omega

theorem nLeft_zero (x lo L m : Nat) (ps : List Part) (h : fits m lo L ps = true) (hx : x ≤ lo) :
    nLeft x ps = 0 := by
  induction ps generalizing lo with
  | nil => rfl
  | cons p t ih =>
    simp [fits] at h
    obtain ⟨⟨h1, _⟩, h3⟩ := h
    have := ih _ h3 (by omega)
    simp [nLeft, this]
    omega

/-- Parts of length ≥ 5 starting before x number at most (x − lo + 4) / 5. -/
theorem nLeft_bound (x lo L : Nat) (ps : List Part) (h : fits 5 lo L ps = true) :
    nLeft x ps = 0 ∨ 5 * nLeft x ps + lo ≤ x + 4 := by
  induction ps generalizing lo with
  | nil => simp [nLeft]
  | cons p t ih =>
    simp [fits] at h
    obtain ⟨⟨h1, h2⟩, h3⟩ := h
    by_cases hp : p.off < x
    · right
      rcases ih _ h3 with h0 | hb
      · simp [nLeft, hp, h0]; omega
      · simp [nLeft, hp]; omega
    · left
      have := nLeft_zero x (p.off + p.len) L 5 t h3 (by omega)
      simp [nLeft, hp, this]

theorem nRight_le (y : Nat) (ps : List Part) : nRight y ps ≤ ps.length := by
  induction ps with
  | nil => simp [nRight]
  | cons p t ih => simp only [nRight, List.length]; split <;> omega

/-- Parts of length ≥ 5 ending after y number at most (L − y + 4) / 5. -/
theorem nRight_bound (y lo L : Nat) (ps : List Part) (h : fits 5 lo L ps = true) :
    nRight y ps = 0 ∨ 5 * nRight y ps + y ≤ L + 4 := by
  induction ps generalizing lo with
  | nil => simp [nRight]
  | cons p t ih =>
    simp [fits] at h
    obtain ⟨⟨h1, h2⟩, h3⟩ := h
    by_cases hp : y < p.off + p.len
    · right
      have hl := fits_len 5 _ L t h3
      have hr := nRight_le y t
      simp [nRight, hp]
      omega
    · rcases ih _ h3 with h0 | hb
      · left; simp [nRight, hp, h0]
      · right; simp [nRight, hp]; omega

theorem kIn_full (lo L m : Nat) (ps : List Part) (h : fits m lo L ps = true) : kIn 0 L ps = kAll ps := by
  induction ps generalizing lo with
  | nil => rfl
  | cons p t ih =>
    simp [fits] at h
    obtain ⟨⟨h1, h2⟩, h3⟩ := h
    have hl := fits_len m _ L t h3
    have := ih _ h3
    simp only [kIn, kAll, this]
    by_cases hc : p.counted = true <;> simp [hc] <;> omega

/-- **Lemma S.** For disjoint parts of length ≥ 5 in [0, L), every span's
bound is at most the full span's. -/
theorem lemmaS (L e : Nat) (ps : List Part) (h : fits 5 0 L ps = true) (x y : Nat) (hy : y ≤ L) :
    spanValue L ps e x y ≤ spanValue L ps e 0 L := by
  have hs := counted_split x y ps
  have hf := kIn_full 0 L 5 ps h
  have hl := nLeft_bound x 0 L ps h
  have hr := nRight_bound y 0 L ps h
  simp only [spanValue, hf]
  by_cases hx : 0 < x <;> by_cases hyL : y < L <;> simp [hx, hyL] <;> omega

/-- Parts of length 4 do not suffice: three counted parts of 4 tiling
[0, 12); the span [9, 12) beats the full span. -/
theorem lemmaS_length_four :
    ∃ L ps x y, fits 4 0 L ps = true ∧ x ≤ y ∧ y ≤ L ∧ spanValue L ps 0 0 L < spanValue L ps 0 x y := by
  refine ⟨12, [⟨0, 4, true⟩, ⟨4, 4, true⟩, ⟨8, 4, true⟩], 9, 12, ?_, ?_, ?_, ?_⟩
  · decide
  · decide
  · decide
  · decide

end Certa.SpanS
