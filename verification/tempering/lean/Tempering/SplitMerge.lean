import Mathlib
import Tempering.SwapKernel

/-!
# Metropolis-Hastings with a general positive proposal, and the sequential-proposal ratio

The split-merge move proposes a whole allocation of the free items by a sequential
(item-by-item) draw and accepts with the Metropolis-Hastings ratio. Two facts are proved here
on a finite state space:

* `mhK_invariant`: the MH kernel with any strictly positive proposal `Q` leaves a strictly
  positive `π` invariant (detailed balance, rows sum to one);
* `seq_ratio`: if the target factorises along the visiting order as `∏ s_t` and the proposal is
  the sequential draw `∏ s_t / Z_t` (each step normalised), the MH ratio is
  `∏ Z_t(proposed) / ∏ Z_t(current)`, the form coded in `mixtureModel.cpp`.

Finite state spaces and exact real arithmetic only; nothing about ergodicity or mixing.
-/

open Finset

set_option linter.unusedSectionVars false

namespace Tempering

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω]

/-- MH acceptance probability for proposing `y` from `x`. -/
noncomputable def mhAcc (π : Ω → ℝ) (Q : Ω → Ω → ℝ) (x y : Ω) : ℝ :=
  min 1 (π y * Q y x / (π x * Q x y))

/-- The MH kernel: proposals accepted with `mhAcc`, rejection mass stays put. -/
noncomputable def mhK (π : Ω → ℝ) (Q : Ω → Ω → ℝ) (x y : Ω) : ℝ :=
  Q x y * mhAcc π Q x y + (if y = x then 1 - ∑ z, Q x z * mhAcc π Q x z else 0)

lemma mhK_row_sum (π : Ω → ℝ) (Q : Ω → Ω → ℝ) (x : Ω) : ∑ y, mhK π Q x y = 1 := by
  unfold mhK
  rw [Finset.sum_add_distrib]
  simp [Finset.sum_ite_eq']

/-- `π x * Q x y * acc x y = min (π x * Q x y) (π y * Q y x)`. -/
lemma flow_eq_min {π : Ω → ℝ} {Q : Ω → Ω → ℝ} (hπ : ∀ x, 0 < π x) (hQ : ∀ x y, 0 < Q x y)
    (x y : Ω) : π x * (Q x y * mhAcc π Q x y) = min (π x * Q x y) (π y * Q y x) := by
  have ha : 0 < π x * Q x y := mul_pos (hπ x) (hQ x y)
  have hb : 0 < π y * Q y x := mul_pos (hπ y) (hQ y x)
  have := mul_min_one_div ha hb
  unfold mhAcc
  calc π x * (Q x y * min 1 (π y * Q y x / (π x * Q x y)))
      = (π x * Q x y) * min 1 (π y * Q y x / (π x * Q x y)) := by ring
    _ = _ := this

theorem mhK_detailedBalance {π : Ω → ℝ} {Q : Ω → Ω → ℝ} (hπ : ∀ x, 0 < π x)
    (hQ : ∀ x y, 0 < Q x y) (x y : Ω) : π x * mhK π Q x y = π y * mhK π Q y x := by
  by_cases hxy : y = x
  · subst hxy; rfl
  · have hxy' : ¬ x = y := fun e => hxy e.symm
    unfold mhK
    rw [if_neg hxy, if_neg hxy', add_zero, add_zero, flow_eq_min hπ hQ, flow_eq_min hπ hQ, min_comm]

/-- The Metropolis-Hastings kernel with a strictly positive proposal leaves `π` invariant. -/
theorem mhK_invariant {π : Ω → ℝ} {Q : Ω → Ω → ℝ} (hπ : ∀ x, 0 < π x) (hQ : ∀ x y, 0 < Q x y) :
    Invariant π (mhK π Q) :=
  invariant_of_detailedBalance (mhK_row_sum π Q) (mhK_detailedBalance hπ hQ)

/-- Sequential proposals: target `∏ s`, proposal `∏ s / Z` (positive factors). The MH ratio
`π(y) q(x) / (π(x) q(y))` is `∏ Z(y) / ∏ Z(x)`. Here `sx, Zx` are the step scores and
normalisers along the path of the current state, `sy, Zy` along the proposed one. -/
theorem seq_ratio {ι : Type*} (T : Finset ι) (sx sy Zx Zy : ι → ℝ)
    (hsx : ∀ t, 0 < sx t) (hsy : ∀ t, 0 < sy t) (hZx : ∀ t, 0 < Zx t) (hZy : ∀ t, 0 < Zy t) :
    ((∏ t ∈ T, sy t) * (∏ t ∈ T, sx t / Zx t)) / ((∏ t ∈ T, sx t) * (∏ t ∈ T, sy t / Zy t))
      = (∏ t ∈ T, Zy t) / (∏ t ∈ T, Zx t) := by
  have e1 : ∏ t ∈ T, sx t / Zx t = (∏ t ∈ T, sx t) / ∏ t ∈ T, Zx t := Finset.prod_div_distrib (f := sx) (g := Zx)
  have e2 : ∏ t ∈ T, sy t / Zy t = (∏ t ∈ T, sy t) / ∏ t ∈ T, Zy t := Finset.prod_div_distrib (f := sy) (g := Zy)
  have p1 : 0 < ∏ t ∈ T, sx t := Finset.prod_pos (fun t _ => hsx t)
  have p2 : 0 < ∏ t ∈ T, sy t := Finset.prod_pos (fun t _ => hsy t)
  have p3 : 0 < ∏ t ∈ T, Zx t := Finset.prod_pos (fun t _ => hZx t)
  have p4 : 0 < ∏ t ∈ T, Zy t := Finset.prod_pos (fun t _ => hZy t)
  rw [e1, e2]
  field_simp

end Tempering
