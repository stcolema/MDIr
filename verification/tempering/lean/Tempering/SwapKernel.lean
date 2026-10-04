import Mathlib

/-!
# Metropolis kernels with a deterministic involutive proposal, on a finite space

This is the finite-state skeleton of the replica-exchange move: a deterministic
proposal `σ` that is an involution (`σ (σ x) = x`) is accepted with probability
`min 1 (π (σ x) / π x)`. We prove that the resulting kernel satisfies detailed
balance with respect to any strictly positive `π`, hence leaves `π` invariant,
and that invariance is preserved by composition and mixing of kernels.

What this does and does not cover: finite state spaces and exact real
arithmetic only. It says nothing about measure-theoretic generalisations,
ergodicity or rates of convergence.
-/

open Finset

set_option linter.unusedSectionVars false
set_option linter.unusedSimpArgs false

namespace Tempering

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω]

/-- `a * min 1 (b / a) = min a b` for positive `a`, `b`. -/
lemma mul_min_one_div {a b : ℝ} (ha : 0 < a) (hb : 0 < b) : a * min 1 (b / a) = min a b := by
  rcases le_total a b with h | h
  · have h1 : 1 ≤ b / a := by rw [le_div_iff₀ ha]; linarith
    rw [min_eq_left h1, min_eq_left h]; ring
  · have h1 : b / a ≤ 1 := by rw [div_le_iff₀ ha]; linarith
    rw [min_eq_right h1, min_eq_right h]; field_simp

/-- Acceptance probability of the deterministic proposal `σ x` from `x`. -/
noncomputable def acc (π : Ω → ℝ) (σ : Ω → Ω) (x : Ω) : ℝ := min 1 (π (σ x) / π x)

lemma acc_nonneg {π : Ω → ℝ} (hπ : ∀ x, 0 < π x) (σ : Ω → Ω) (x : Ω) : 0 ≤ acc π σ x :=
  le_min zero_le_one (div_nonneg (hπ _).le (hπ x).le)

lemma acc_le_one (π : Ω → ℝ) (σ : Ω → Ω) (x : Ω) : acc π σ x ≤ 1 := min_le_left _ _

/-- The Metropolis kernel: move to `σ x` with probability `acc`, otherwise stay. -/
noncomputable def swapK (π : Ω → ℝ) (σ : Ω → Ω) (x y : Ω) : ℝ :=
  (if y = σ x then acc π σ x else 0) + (if y = x then 1 - acc π σ x else 0)

lemma swapK_nonneg {π : Ω → ℝ} (hπ : ∀ x, 0 < π x) (σ : Ω → Ω) (x y : Ω) :
    0 ≤ swapK π σ x y := by
  unfold swapK
  have h1 := acc_nonneg hπ σ x
  have h2 := acc_le_one π σ x
  split_ifs <;> linarith

/-- Rows sum to one: it is a Markov kernel. -/
lemma swapK_row_sum (π : Ω → ℝ) (σ : Ω → Ω) (x : Ω) : ∑ y, swapK π σ x y = 1 := by
  unfold swapK
  rw [Finset.sum_add_distrib]
  simp [Finset.sum_ite_eq']

/-- Detailed balance of the Metropolis kernel with an involutive proposal. -/
theorem swapK_detailedBalance (π : Ω → ℝ) (hπ : ∀ x, 0 < π x) (σ : Ω → Ω)
    (hσ : ∀ x, σ (σ x) = x) (x y : Ω) :
    π x * swapK π σ x y = π y * swapK π σ y x := by
  by_cases hxy : y = x
  · subst hxy; rfl
  · have hiff : y = σ x ↔ x = σ y := by
      constructor
      · intro h; rw [h, hσ]
      · intro h; rw [h, hσ]
    by_cases h : y = σ x
    · have h' : x = σ y := hiff.mp h
      have hxy' : ¬ x = y := fun e => hxy e.symm
      have e1 : π x * swapK π σ x y = min (π x) (π y) := by
        unfold swapK
        rw [if_pos h, if_neg hxy, add_zero]
        unfold acc
        rw [← h]
        exact mul_min_one_div (hπ x) (hπ y)
      have e2 : π y * swapK π σ y x = min (π y) (π x) := by
        unfold swapK
        rw [if_pos h', if_neg hxy', add_zero]
        unfold acc
        rw [← h']
        exact mul_min_one_div (hπ y) (hπ x)
      rw [e1, e2, min_comm]
    · have h' : ¬ x = σ y := fun e => h (hiff.mpr e)
      have hxy' : ¬ x = y := fun e => hxy e.symm
      unfold swapK
      rw [if_neg h, if_neg hxy, if_neg h', if_neg hxy']
      simp

/-- `π` is invariant for the kernel `K` (a matrix indexed by states). -/
def Invariant (π : Ω → ℝ) (K : Ω → Ω → ℝ) : Prop := ∀ y, ∑ x, π x * K x y = π y

/-- Detailed balance with rows summing to one gives invariance. -/
theorem invariant_of_detailedBalance {π : Ω → ℝ} {K : Ω → Ω → ℝ}
    (hrow : ∀ x, ∑ y, K x y = 1) (hdb : ∀ x y, π x * K x y = π y * K y x) : Invariant π K := by
  intro y
  calc ∑ x, π x * K x y = ∑ x, π y * K y x := Finset.sum_congr rfl (fun x _ => hdb x y)
    _ = π y * ∑ x, K y x := by rw [Finset.mul_sum]
    _ = π y := by rw [hrow y, mul_one]

/-- The swap kernel leaves `π` invariant. -/
theorem swapK_invariant (π : Ω → ℝ) (hπ : ∀ x, 0 < π x) (σ : Ω → Ω)
    (hσ : ∀ x, σ (σ x) = x) : Invariant π (swapK π σ) :=
  invariant_of_detailedBalance (swapK_row_sum π σ) (swapK_detailedBalance π hπ σ hσ)

/-- Composition of kernels (apply `K₁`, then `K₂`). -/
def comp (K₁ K₂ : Ω → Ω → ℝ) (x z : Ω) : ℝ := ∑ y, K₁ x y * K₂ y z

/-- Invariance is preserved by composition: so any fixed or alternating
sequence of invariant kernels (as in even/odd exchange rounds) is invariant. -/
theorem invariant_comp {π : Ω → ℝ} {K₁ K₂ : Ω → Ω → ℝ}
    (h₁ : Invariant π K₁) (h₂ : Invariant π K₂) : Invariant π (comp K₁ K₂) := by
  intro z
  calc ∑ x, π x * comp K₁ K₂ x z
      = ∑ x, ∑ y, π x * K₁ x y * K₂ y z := by
        refine Finset.sum_congr rfl (fun x _ => ?_)
        unfold comp; rw [Finset.mul_sum]
        exact Finset.sum_congr rfl (fun y _ => by ring)
    _ = ∑ y, ∑ x, π x * K₁ x y * K₂ y z := Finset.sum_comm
    _ = ∑ y, (∑ x, π x * K₁ x y) * K₂ y z := by
        refine Finset.sum_congr rfl (fun y _ => ?_)
        rw [Finset.sum_mul]
    _ = ∑ y, π y * K₂ y z := by
        refine Finset.sum_congr rfl (fun y _ => ?_)
        rw [h₁ y]
    _ = π z := h₂ z

/-- A convex mixture of invariant kernels is invariant (stochastic even-odd
exchange chooses the parity with probability one half). -/
theorem invariant_mix {π : Ω → ℝ} {K₁ K₂ : Ω → Ω → ℝ} (p : ℝ)
    (h₁ : Invariant π K₁) (h₂ : Invariant π K₂) : Invariant π (fun x y => p * K₁ x y + (1 - p) * K₂ x y) := by
  intro y
  calc ∑ x, π x * (p * K₁ x y + (1 - p) * K₂ x y)
      = p * ∑ x, π x * K₁ x y + (1 - p) * ∑ x, π x * K₂ x y := by
        rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
        exact Finset.sum_congr rfl (fun x _ => by ring)
    _ = π y := by rw [h₁ y, h₂ y]; ring

/-- The identity kernel is invariant (the base case for an empty round). -/
theorem invariant_id (π : Ω → ℝ) : Invariant π (fun x y => if y = x then 1 else 0) := by
  intro y
  simp [Finset.sum_ite_eq']

end Tempering
