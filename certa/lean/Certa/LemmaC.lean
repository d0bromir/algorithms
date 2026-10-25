import Certa.LemmaW

/-
Lemma C (chains; THEORY.md §4), its geometric core. Occurrences of the parts
are sorted by diagonal (reference position minus read position) and split
into chains where consecutive diagonals differ by more than pad. For an
alignment A with gap total G ≤ pad:
- `diag_drift`: two aligned bases of A lie on diagonals at most G apart (a
  diagonal moves only at gaps: −1 per inserted base, +k per deletion of k);
  with Lemma W, at most pad(b) apart when A scores above b (`diag_within_pad`);
- `no_split_within`: no split falls between two diagonals at most pad
  apart; so A's exact parts' occurrences lie in one chain, and A's aligned
  bases on that chain's diagonals ± pad.
-/
namespace Certa

def stIns : St → Int
  | .i => 1
  | _ => 0

/-- How far the diagonal moves over the bases: +deleted − inserted. -/
def drift : List Base → Int
  | [] => 0
  | b :: t => (b.del : Int) - stIns b.s + drift t

theorem drift_append (u v : List Base) : drift (u ++ v) = drift u + drift v := by
  induction u with
  | nil => simp [drift]
  | cons b t ih => simp [drift, ih]; omega

theorem gaps_append (u v : List Base) : gaps (u ++ v) = gaps u + gaps v := by
  induction u with
  | nil => simp [gaps]
  | cons b t ih => simp [gaps, ih]; omega

theorem drift_le_gaps : ∀ l : List Base, drift l ≤ gaps l ∧ -(gaps l : Int) ≤ drift l
  | [] => by simp [drift, gaps]
  | b :: t => by
    have := drift_le_gaps t
    obtain ⟨s, d⟩ := b
    cases s <;> simp [drift, gaps, stIns, stGap] <;> omega

/-- **Diagonal drift.** The bases after a prefix u and after a longer prefix
u ++ v of an alignment's transcript lie on diagonals that differ by
drift v, at most the gap total. -/
theorem diag_drift (u v w : List Base) :
    (drift (u ++ v) - drift u) ≤ gaps (u ++ v ++ w) ∧
      -(gaps (u ++ v ++ w) : Int) ≤ drift (u ++ v) - drift u := by
  have h1 := drift_append u v
  have h2 := drift_le_gaps v
  have h3 := gaps_append (u ++ v) w
  have h4 := gaps_append u v
  omega

/-- With Lemma W: in an alignment scoring above b, any two aligned bases lie
on diagonals at most pad(b) = max(0, L − b − 7) apart. -/
theorem diag_within_pad (c₁ c₂ : Nat) (u v w : List Base) (b : Int)
    (hs : localScore c₁ c₂ (u ++ v ++ w) > b) :
    (drift (u ++ v) - drift u) ≤ max 0 ((c₁ + c₂ + span (u ++ v ++ w) : Nat) - b - 7) ∧
      -(max 0 ((c₁ + c₂ + span (u ++ v ++ w) : Nat) - b - 7)) ≤ drift (u ++ v) - drift u := by
  have hw := lemmaW c₁ c₂ (u ++ v ++ w) b hs
  have hd := diag_drift u v w
  omega

/-- **No split between close diagonals.** In sorted diagonals ds, a split
between positions k − 1 and k (ds k − ds (k − 1) > pad) with i < k ≤ j makes
ds j − ds i > pad: two occurrences at most pad apart are in one chain. -/
theorem no_split_within (ds : Nat → Int) (mono : ∀ a b, a ≤ b → ds a ≤ ds b) (pad : Int)
    (i j k : Nat) (hik : i < k) (hkj : k ≤ j) (hsplit : ds k - ds (k - 1) > pad) :
    ds j - ds i > pad := by
  have h1 := mono k j hkj
  have h2 := mono i (k - 1) (by omega)
  omega

/-- Parts marked free (the parts of chain C) that cover every exactly
occurring part leave at least as many inexact parts as there are parts not
free. -/
theorem notFree_le_inexact :
    ∀ segs : List (Bool × (List Base × List Base)),
      (∀ x ∈ segs, exact x.2.2 = true → x.1 = true) →
      (segs.filter (fun x => !x.1)).length ≤
        ((segs.map Prod.snd).filter (fun gp => !exact gp.2)).length
  | [], _ => by simp
  | x :: t, h => by
    have ih := notFree_le_inexact t (fun y hy => h y (by simp [hy]))
    obtain ⟨f, g, p⟩ := x
    have hx : exact p = true → f = true := h (f, g, p) (by simp)
    cases f <;> cases hp : exact p <;> simp_all <;> omega

/-- **Lemma C's bound** (σ(A) ≤ ub(C)). If every exactly occurring part of
an alignment is among the free parts (the parts of the chain C it lies in),
it loses at least 5 for each part of its span that is not free: Lemma L with
C's parts treated as free. -/
theorem lemmaC_bound (segs : List (Bool × (List Base × List Base))) (tail : List Base)
    (hlen : ∀ x ∈ segs, 3 ≤ x.2.2.length)
    (hfree : ∀ x ∈ segs, exact x.2.2 = true → x.1 = true) :
    5 * (segs.filter (fun x => !x.1)).length ≤ loss (assemble (segs.map Prod.snd) tail) := by
  have hL := lemmaL (segs.map Prod.snd) tail (by
    intro gp hgp
    obtain ⟨x, hx, rfl⟩ := List.mem_map.mp hgp
    exact hlen x hx)
  have hc := notFree_le_inexact segs hfree
  omega

end Certa
