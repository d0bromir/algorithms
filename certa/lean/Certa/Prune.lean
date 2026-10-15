/-
Pruning the band (THEORY.md "Pruning the band"; certa/longread.h band_best
with a floor). The banded DP of the long-read certificates needs only
alignments above a floor. A cell of value v at read base f leads at most to
v + (L − 1 − f), and is dead when that is ≤ the floor. Two facts make the
pruned DP exact above the floor; here they are proved for alignment paths
(the DP's cell value is the best path ending there):
- `prefix_live`: every cell on a path to a live cell is live, so dropping
  dead cells changes no live cell's value;
- `above_live`: once no new piece can lead above the floor, the cell above a
  live cell (one read base earlier, same diagonal) is live, so each row's
  live cells lie in the previous row's live range.

A path: a first aligned pair (read base f₀, with start value g f₀ = 0 for
f₀ = 0, else −5 for the clipped start), then aligned pairs (+1 a match, −4 a
mismatch), inserted read bases and deleted reference bases, a gap run of k
costing 6 + k.
-/
namespace Certa.Prune

inductive Op where
  | p (hit : Bool)  -- an aligned pair: a match or a mismatch
  | i               -- an inserted read base
  | d               -- a deleted reference base
  deriving DecidableEq, Repr

/-- The kind of the previous op: a gap of the same kind extends (−1), any
other gap opens (−7). -/
inductive Kind where
  | pair | ins | del
  deriving DecidableEq, Repr

def pairScore (hit : Bool) : Int := if hit then 1 else -4

def sc : Kind → List Op → Int
  | _, [] => 0
  | _, .p b :: t => pairScore b + sc .pair t
  | k, .i :: t => (if k = .ins then -1 else -7) + sc .ins t
  | k, .d :: t => (if k = .del then -1 else -7) + sc .del t

def lastKind : Kind → List Op → Kind
  | k, [] => k
  | _, .p _ :: t => lastKind .pair t
  | _, .i :: t => lastKind .ins t
  | _, .d :: t => lastKind .del t

/-- Read bases (pairs and inserted bases) and reference bases (pairs and
deleted bases) of the ops. -/
def reads : List Op → Nat
  | [] => 0
  | .p _ :: t => 1 + reads t
  | .i :: t => 1 + reads t
  | .d :: t => reads t

def refs : List Op → Nat
  | [] => 0
  | .p _ :: t => 1 + refs t
  | .i :: t => refs t
  | .d :: t => 1 + refs t

def pairs : List Op → Nat
  | [] => 0
  | .p _ :: t => 1 + pairs t
  | _ :: t => pairs t

theorem sc_append (k : Kind) (u v : List Op) : sc k (u ++ v) = sc k u + sc (lastKind k u) v := by
  induction u generalizing k with
  | nil => simp [sc, lastKind]
  | cons o t ih => cases o <;> simp [sc, lastKind, ih] <;> omega

/-- After a pair, every gap opens: no previous kind does worse. -/
theorem sc_pair_le (k : Kind) : ∀ l, sc .pair l ≤ sc k l
  | [] => Int.le_refl _
  | .p _ :: _ => Int.le_refl _
  | .i :: _ => by cases k <;> simp [sc]
  | .d :: _ => by cases k <;> simp [sc]

/-- Each read base adds at most 1; a deleted base adds nothing. -/
theorem sc_le_reads (k : Kind) : ∀ l, sc k l ≤ reads l
  | [] => Int.le_refl _
  | .p b :: t => by
    have := sc_le_reads .pair t
    simp [sc, reads, pairScore]; split <;> omega
  | .i :: t => by
    have := sc_le_reads .ins t
    simp [sc, reads]; split <;> omega
  | .d :: t => by
    have := sc_le_reads .del t
    simp [sc, reads]; split <;> omega

/-- A path: its start value, its first pair, then the rest. -/
structure Path where
  g : Int
  hit : Bool
  rest : List Op

def Path.score (P : Path) : Int := P.g + pairScore P.hit + sc .pair P.rest
/-- Read bases after the first; the path ends on row f₀ + P.len. -/
def Path.len (P : Path) : Nat := reads P.rest
/-- Its last cell's diagonal, relative to its first. -/
def Path.diag (P : Path) : Int := (refs P.rest : Int) - reads P.rest

theorem reads_append (u v : List Op) : reads (u ++ v) = reads u + reads v := by
  induction u with
  | nil => simp [reads]
  | cons o t ih => cases o <;> simp [reads, ih] <;> omega

theorem refs_append (u v : List Op) : refs (u ++ v) = refs u + refs v := by
  induction u with
  | nil => simp [refs]
  | cons o t ih => cases o <;> simp [refs, ih] <;> omega

/-- **Prefixes.** The path through its first ops u has a value and a bound
at least the whole path's: a cell on a path to a live cell is live. -/
theorem prefix_live (g : Int) (hit : Bool) (u v : List Op) (L f₀ : Nat) (floor : Int)
    (hlive : g + pairScore hit + sc .pair (u ++ v) + ((L : Int) - 1 - (f₀ + reads (u ++ v))) > floor) :
    g + pairScore hit + sc .pair u + ((L : Int) - 1 - (f₀ + reads u)) > floor := by
  have hreads := reads_append u v
  have h1 := sc_append .pair u v
  have h2 := sc_le_reads (lastKind .pair u) v
  omega

/-- **Removing a pair.** Dropping an aligned pair after the first keeps the
diagonal, ends one read base earlier, and loses at most 1. -/
theorem remove_pair (P : Path) (u v : List Op) (b : Bool) (h : P.rest = u ++ .p b :: v) :
    let P' : Path := ⟨P.g, P.hit, u ++ v⟩
    P'.score ≥ P.score - 1 ∧ P'.len + 1 = P.len ∧ P'.diag = P.diag := by
  intro P'
  have e1 := sc_append .pair u (.p b :: v)
  have e2 := sc_append .pair u v
  have hle := sc_pair_le (lastKind .pair u) v
  have hr := reads_append u (.p b :: v)
  have hr' := reads_append u v
  have hf := refs_append u (.p b :: v)
  have hf' := refs_append u v
  simp only [Path.score, Path.len, Path.diag, P', h] at *
  simp only [sc, reads, refs, pairScore] at e1 hr hf
  refine ⟨?_, ?_, ?_⟩
  · split at e1 <;> omega
  · omega
  · omega

/-- The start value of a piece starting at read base f₀. -/
def g (f₀ : Nat) : Int := if f₀ = 0 then 0 else -5

/-- A path of one pair and then gaps only. -/
def gapsOnly : List Op → Prop
  | [] => True
  | .p _ :: _ => False
  | _ :: t => gapsOnly t

theorem gaps_bound (k : Kind) : ∀ l, gapsOnly l → l ≠ [] → k = .pair → sc k l + 6 + reads l ≤ 0
  | [], _, h, _ => absurd rfl h
  | .p _ :: _, hg, _, _ => absurd hg (by simp [gapsOnly])
  | .i :: t, hg, _, hk => by
    subst hk
    have := sc_le_reads .ins t
    have hgt : gapsOnly t := hg
    -- t's ops are gaps: each costs at least 1, so sc .ins t ≤ -(reads t)
    have hneg := gaps_neg .ins t hgt
    simp [sc, reads]; omega
  | .d :: t, hg, _, hk => by
    subst hk
    have hgt : gapsOnly t := hg
    have hneg := gaps_neg .del t hgt
    simp [sc, reads]; omega
where
  /-- Gaps cost at least 1 per base: sc ≤ −(reads). -/
  gaps_neg (k : Kind) : ∀ l, gapsOnly l → sc k l + reads l ≤ 0
    | [], _ => by simp [sc, reads]
    | .p _ :: _, hg => absurd hg (by simp [gapsOnly])
    | .i :: t, hg => by
      have := gaps_neg .ins t hg
      simp [sc, reads]; split <;> omega
    | .d :: t, hg => by
      have := gaps_neg .del t hg
      simp [sc, reads]; split <;> omega

/-- A path with one aligned pair, or more. -/
theorem one_or_two (rest : List Op) : gapsOnly rest ∨ ∃ u b v, rest = u ++ .p b :: v := by
  induction rest with
  | nil => exact Or.inl trivial
  | cons o t ih =>
    cases o with
    | p b => exact Or.inr ⟨[], b, t, rfl⟩
    | i =>
      rcases ih with h | ⟨u, b, v, h⟩
      · exact Or.inl h
      · exact Or.inr ⟨.i :: u, b, v, by simp [h]⟩
    | d =>
      rcases ih with h | ⟨u, b, v, h⟩
      · exact Or.inl h
      · exact Or.inr ⟨.d :: u, b, v, by simp [h]⟩

/-- **The cell above.** A path from read base f₀ to row f = f₀ + len, live
(its value plus the L − 1 − f bases after f exceeds the floor), in a row
where no new piece can lead above the floor (g f + 1 + L − 1 − f ≤ floor):
it has two aligned pairs, and dropping one gives a path from f₀ to row f − 1
on the same diagonal, also live. So row f's live cells lie in the previous
row's live range. -/
theorem above_live (P : Path) (f₀ L : Nat) (floor : Int) (hg : P.g = g f₀)
    (hlive : P.score + ((L : Int) - 1 - (f₀ + P.len)) > floor)
    (hnew : g (f₀ + P.len) + 1 + ((L : Int) - 1 - (f₀ + P.len)) ≤ floor) :
    ∃ P' : Path, P'.g = P.g ∧ P'.len + 1 = P.len ∧ P'.diag = P.diag ∧
      P'.score + ((L : Int) - 1 - (f₀ + P'.len)) > floor := by
  rcases one_or_two P.rest with hgo | ⟨u, b, v, h⟩
  · -- one pair and gaps: already below a new piece at row f
    exfalso
    have hps : pairScore P.hit ≤ 1 := by simp [pairScore]; split <;> omega
    have hg0 : g f₀ ≤ 0 := by simp [g]; split <;> omega
    by_cases hnil : P.rest = []
    · have hl : P.len = 0 := by simp [Path.len, hnil, reads]
      have hs : P.score ≤ g f₀ + 1 := by
        simp only [Path.score, hnil, sc]; omega
      rw [hl] at hlive hnew
      simp only [Nat.add_zero] at hlive hnew
      omega
    · have hb := gaps_bound .pair P.rest hgo hnil rfl
      have hgf : g (f₀ + P.len) ≥ -5 := by simp [g]; split <;> omega
      simp only [Path.score, Path.len] at hlive hnew hgf
      omega
  · obtain ⟨h1, h2, h3⟩ := remove_pair P u v b h
    refine ⟨⟨P.g, P.hit, u ++ v⟩, rfl, h2, h3, ?_⟩
    have : ((⟨P.g, P.hit, u ++ v⟩ : Path).len : Int) + 1 = P.len := by exact_mod_cast h2
    omega

end Certa.Prune
