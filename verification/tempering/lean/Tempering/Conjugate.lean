import Mathlib
import Tempering.SwapKernel

/-!
# Tempered full conditionals

1. Heat-bath (Gibbs) kernels on a finite state space: redrawing a block from
   its full conditional leaves the target invariant. Applied to the tempered
   target `π_β ∝ L^β P`, this is why a sampler that draws each block from its
   *tempered* full conditional targets `π_β`.
2. The completing-the-square identity behind the tempered Normal-inverse-gamma
   (and Normal-inverse-Wishart) update: the likelihood `L^β` of `n`
   observations updates the prior exactly as `β n` observations with `β` times
   the centred sum of squares.
3. The tempered categorical likelihood adds `β n_k` to a Dirichlet exponent.

Scalar versions only for 2 and 3: the matrix identity for the
normal-inverse-Wishart case is checked symbolically with SymPy and by
simulation (see `verification/tempering`).
-/

open Finset

set_option linter.unusedSectionVars false
set_option linter.unusedSimpArgs false

namespace Tempering

section Gibbs

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω] (r : Ω → Ω → Prop) [DecidableRel r]

/-- Total target mass of the block of `x` (the block is the equivalence class). -/
noncomputable def blockMass (π : Ω → ℝ) (x : Ω) : ℝ :=
  ∑ z ∈ Finset.univ.filter (fun z => r z x), π z

/-- Heat-bath kernel: redraw the state within its block with probability
proportional to the target. -/
noncomputable def gibbsK (π : Ω → ℝ) (x y : Ω) : ℝ :=
  if r y x then π y / blockMass r π x else 0

variable {r}

lemma blockMass_eq (hsymm : ∀ x y, r x y → r y x) (htrans : ∀ x y z, r x y → r y z → r x z)
    {π : Ω → ℝ} {x y : Ω} (h : r y x) : blockMass r π x = blockMass r π y := by
  unfold blockMass
  refine Finset.sum_congr ?_ (fun _ _ => rfl)
  ext z
  simp only [Finset.mem_filter, Finset.mem_univ, true_and]
  constructor
  · intro hz; exact htrans z x y hz (hsymm y x h)
  · intro hz; exact htrans z y x hz h

lemma blockMass_pos (hrefl : ∀ x, r x x) {π : Ω → ℝ} (hπ : ∀ x, 0 < π x) (x : Ω) :
    0 < blockMass r π x := by
  unfold blockMass
  apply Finset.sum_pos (fun z _ => hπ z)
  exact ⟨x, by simp [hrefl x]⟩

lemma gibbsK_row_sum (hrefl : ∀ x, r x x) {π : Ω → ℝ} (hπ : ∀ x, 0 < π x) (x : Ω) :
    ∑ y, gibbsK r π x y = 1 := by
  unfold gibbsK
  have hm := blockMass_pos hrefl hπ x
  rw [← Finset.sum_filter, ← Finset.sum_div]
  exact div_self hm.ne'

theorem gibbsK_detailedBalance (hsymm : ∀ x y, r x y → r y x)
    (htrans : ∀ x y z, r x y → r y z → r x z) {π : Ω → ℝ} (x y : Ω) :
    π x * gibbsK r π x y = π y * gibbsK r π y x := by
  unfold gibbsK
  by_cases h : r y x
  · have h' : r x y := hsymm y x h
    rw [if_pos h, if_pos h', blockMass_eq hsymm htrans h]
    ring
  · have h' : ¬ r x y := fun e => h (hsymm x y e)
    rw [if_neg h, if_neg h']; ring

/-- Drawing a block from its full conditional leaves the target invariant. -/
theorem gibbsK_invariant (hrefl : ∀ x, r x x) (hsymm : ∀ x y, r x y → r y x)
    (htrans : ∀ x y z, r x y → r y z → r x z) {π : Ω → ℝ} (hπ : ∀ x, 0 < π x) :
    Invariant π (gibbsK r π) :=
  invariant_of_detailedBalance (gibbsK_row_sum hrefl hπ) (gibbsK_detailedBalance hsymm htrans)

end Gibbs

/-- Tempered completing the square. With `n` observations of sum `S1` and sum of
squares `S2`, prior precision weight `κ` and prior mean `ξ`:
`β Σ(xᵢ-μ)² + κ(μ-ξ)²  =  (κ+βn)(μ-μₙ)² + β·SS + (κ βn/(κ+βn))(x̄-ξ)²`
with `μₙ = (κξ + β S1)/(κ+βn)` and `SS = S2 - S1²/n`: the update of a prior by the
likelihood raised to `β` is the untempered update with `βn` observations. -/
theorem tempered_square (β κ ξ μ n S1 S2 : ℝ) (hn : 0 < n) (hk : 0 < κ + β * n) :
    β * (S2 - 2 * μ * S1 + n * μ ^ 2) + κ * (μ - ξ) ^ 2
      = (κ + β * n) * (μ - (κ * ξ + β * S1) / (κ + β * n)) ^ 2
        + β * (S2 - S1 ^ 2 / n)
        + (κ * (β * n) / (κ + β * n)) * (S1 / n - ξ) ^ 2 := by
  field_simp
  ring

/-- The same, stated through the effective sample size `nₑ = β n`: the
posterior mean is the weighted average with `nₑ` and the centred term has the
untempered form with `nₑ`. -/
theorem tempered_square_eff (β κ ξ μ n S1 S2 : ℝ) (hn : 0 < n) (hβ : β ≠ 0) (hk : 0 < κ + β * n) :
    β * (S2 - 2 * μ * S1 + n * μ ^ 2) + κ * (μ - ξ) ^ 2
      = (κ + β * n) * (μ - (κ * ξ + (β * n) * (S1 / n)) / (κ + β * n)) ^ 2
        + β * (S2 - S1 ^ 2 / n)
        + (κ * (β * n) / (κ + β * n)) * (S1 / n - ξ) ^ 2 := by
  have : (β * n) * (S1 / n) = β * S1 := by field_simp
  rw [this]
  exact tempered_square β κ ξ μ n S1 S2 hn hk

/-- Tempered categorical likelihood: `(θ ^ n) ^ β = θ ^ (n β)`, so a Dirichlet
exponent `α - 1` becomes `α - 1 + β n`. -/
theorem tempered_power {θ : ℝ} (hθ : 0 < θ) (n : ℕ) (β : ℝ) :
    (θ ^ n) ^ β = θ ^ ((n : ℝ) * β) := by
  rw [← Real.rpow_natCast, ← Real.rpow_mul hθ.le]

end Tempering
