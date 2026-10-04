import Mathlib
import Tempering.SwapKernel

/-!
# Annealed importance sampling on a finite space

A chain of annealing steps: unnormalised targets `g₀, g₁, …, g_K` (strictly
positive) and kernels `M_k` with `g_k M_k = g_k` (invariance). The sampler draws
`x₀ ∼ g₀ / Z₀`, and for `k = 1..K` multiplies the weight by `g_k(x_{k-1}) / g_{k-1}(x_{k-1})`
and moves `x_{k-1} → x_k` with `M_k`. The unnormalised measure of "weight times
indicator of the current state" evolves by `stepMeasure` (this is the forward
equation of the weighted path; the claim that it is the expectation of
`weight * f(state)` is the definition of the forward measure).

Proved here:

* `annealMeasure_eq_last`: if the start is `g₀ / Z₀` the measure after `K` steps is
  `g_K / Z₀`. Hence `E[W] = Z_K / Z₀` and `E[W f(x_K)] = (1/Z₀) Σ g_K f`:
  the weights are exactly unbiased for the unnormalised target, for ANY number of
  particles and whether or not the kernels mix (finite-sample unbiasedness).
* `bias_bound`: if the start is `μ₀` instead of `π₀ = g₀ / Z₀`, and each incremental
  ratio is at most `R_k`, the bias of `E[W f(x_K)]` for `|f| ≤ B` is at most
  `B * ∏ R_k * Σ |μ₀ - π₀|` (twice the total-variation distance of the start).
  This is the only rigorous statement available for a start that is not exact, and
  `Σ |μ₀ - π₀|` is not estimable from the sampler's output.

Not covered: resampling (checked in exact rational arithmetic for tiny cases,
`sympy_smc_exact.py`, and cited: Del Moral, 2004, Thm 7.4.2), the law of large
numbers and central limit theorem for the self-normalised estimator (cited).
-/

open Finset

set_option linter.unusedSectionVars false
set_option linter.unusedSimpArgs false

namespace Tempering

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω]

/-- One annealing step. -/
structure AnnealStep (Ω : Type*) where
  gprev : Ω → ℝ
  gnext : Ω → ℝ
  M : Ω → Ω → ℝ

/-- Weighted state measure after a step: weight by the incremental ratio, then move. -/
noncomputable def stepMeasure (s : AnnealStep Ω) (ν : Ω → ℝ) (y : Ω) : ℝ :=
  ∑ x, ν x * (s.gnext x / s.gprev x) * s.M x y

/-- The measure after a list of steps. -/
noncomputable def annealMeasure : List (AnnealStep Ω) → (Ω → ℝ) → (Ω → ℝ)
  | [], ν => ν
  | (s :: rest), ν => annealMeasure rest (stepMeasure s ν)

/-- The steps form a valid chain starting from `g`: each step's previous target is
the current one, targets are positive, and each kernel leaves its next target invariant. -/
def ValidChain : (Ω → ℝ) → List (AnnealStep Ω) → Prop
  | _, [] => True
  | g, (s :: rest) =>
      s.gprev = g ∧ (∀ x, 0 < s.gprev x) ∧ Invariant s.gnext s.M ∧ ValidChain s.gnext rest

/-- The final target of a chain. -/
def lastTarget : (Ω → ℝ) → List (AnnealStep Ω) → (Ω → ℝ)
  | g, [] => g
  | _, (s :: rest) => lastTarget s.gnext rest

theorem stepMeasure_eq (s : AnnealStep Ω) (hpos : ∀ x, 0 < s.gprev x)
    (hinv : Invariant s.gnext s.M) : stepMeasure s s.gprev = s.gnext := by
  funext y
  unfold stepMeasure
  calc ∑ x, s.gprev x * (s.gnext x / s.gprev x) * s.M x y = ∑ x, s.gnext x * s.M x y := by
        refine Finset.sum_congr rfl (fun x _ => ?_)
        have := (hpos x).ne'
        field_simp
    _ = s.gnext y := hinv y

/-- Unbiasedness of the weighted measure: started at `g`, it ends at the last target. -/
theorem annealMeasure_eq_last : ∀ (steps : List (AnnealStep Ω)) (g : Ω → ℝ),
    ValidChain g steps → annealMeasure steps g = lastTarget g steps
  | [], g, _ => rfl
  | (s :: rest), g, h => by
      obtain ⟨hg, hpos, hinv, hrest⟩ := h
      have h1 : stepMeasure s g = s.gnext := by
        rw [← hg]; exact stepMeasure_eq s hpos hinv
      show annealMeasure rest (stepMeasure s g) = lastTarget s.gnext rest
      rw [h1]
      exact annealMeasure_eq_last rest s.gnext hrest

theorem stepMeasure_smul (s : AnnealStep Ω) (c : ℝ) (ν : Ω → ℝ) :
    stepMeasure s (fun x => c * ν x) = fun y => c * stepMeasure s ν y := by
  funext y
  unfold stepMeasure
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl (fun x _ => by ring)

theorem annealMeasure_smul : ∀ (steps : List (AnnealStep Ω)) (c : ℝ) (ν : Ω → ℝ),
    annealMeasure steps (fun x => c * ν x) = fun y => c * annealMeasure steps ν y
  | [], c, ν => rfl
  | (s :: rest), c, ν => by
      show annealMeasure rest (stepMeasure s (fun x => c * ν x)) = _
      rw [stepMeasure_smul]
      exact annealMeasure_smul rest c (stepMeasure s ν)

/-- Started from the normalised `g₀ / Z₀`, the weighted measure after `K` steps is
`g_K / Z₀`; in particular `E[W f(x_K)] = (Σ g_K f) / Z₀` for every `f` and every
number of particles, and `E[W] = Z_K / Z₀`. -/
theorem ais_unbiased (steps : List (AnnealStep Ω)) (g : Ω → ℝ) (Z₀ : ℝ) (hZ : Z₀ ≠ 0)
    (h : ValidChain g steps) (f : Ω → ℝ) :
    ∑ y, annealMeasure steps (fun x => g x / Z₀) y * f y = (∑ y, lastTarget g steps y * f y) / Z₀ := by
  have e : (fun x => g x / Z₀) = fun x => (1 / Z₀) * g x := by
    funext x; field_simp
  rw [e, annealMeasure_smul, annealMeasure_eq_last steps g h]
  rw [Finset.sum_div]
  exact Finset.sum_congr rfl (fun y _ => by field_simp)

/-! ### Bias of a start that is not exact -/

/-- The conditional expected weighted test function from state `x` after the steps. -/
noncomputable def annealTest : List (AnnealStep Ω) → (Ω → ℝ) → (Ω → ℝ)
  | [], f => f
  | (s :: rest), f => fun x => (s.gnext x / s.gprev x) * ∑ y, s.M x y * annealTest rest f y

/-- Duality: integrating the weighted measure against `f` equals integrating the
start against the conditional expected weighted test function. -/
theorem annealMeasure_dual : ∀ (steps : List (AnnealStep Ω)) (ν f : Ω → ℝ),
    ∑ y, annealMeasure steps ν y * f y = ∑ x, ν x * annealTest steps f x
  | [], ν, f => rfl
  | (s :: rest), ν, f => by
      show ∑ y, annealMeasure rest (stepMeasure s ν) y * f y = _
      rw [annealMeasure_dual rest (stepMeasure s ν) f]
      have hcons : annealTest (s :: rest) f
          = fun x => (s.gnext x / s.gprev x) * ∑ y, s.M x y * annealTest rest f y := rfl
      rw [hcons]
      unfold stepMeasure
      simp_rw [Finset.sum_mul]
      rw [Finset.sum_comm]
      refine Finset.sum_congr rfl (fun x _ => ?_)
      rw [Finset.mul_sum, Finset.mul_sum]
      refine Finset.sum_congr rfl (fun y _ => ?_)
      ring

/-- Product of the per-step ratio bounds. -/
def ratioBound (R : AnnealStep Ω → ℝ) : List (AnnealStep Ω) → ℝ
  | [] => 1
  | (s :: rest) => R s * ratioBound R rest

theorem annealTest_bound (R : AnnealStep Ω → ℝ) :
    ∀ (steps : List (AnnealStep Ω)) (f : Ω → ℝ) (B : ℝ),
      (∀ s ∈ steps, 0 ≤ R s) →
      (∀ s ∈ steps, ∀ x, 0 ≤ s.gnext x / s.gprev x ∧ s.gnext x / s.gprev x ≤ R s) →
      (∀ s ∈ steps, ∀ x y, 0 ≤ s.M x y) →
      (∀ s ∈ steps, ∀ x, ∑ y, s.M x y = 1) →
      (∀ x, |f x| ≤ B) → ∀ x, |annealTest steps f x| ≤ B * ratioBound R steps
  | [], f, B, _, _, _, _, hf, x => by simpa [annealTest, ratioBound] using hf x
  | (s :: rest), f, B, hR, hr, hM, hrow, hf, x => by
      have hs : s ∈ s :: rest := List.mem_cons_self
      have hrest := annealTest_bound R rest f B
        (fun t ht => hR t (List.mem_cons_of_mem _ ht))
        (fun t ht => hr t (List.mem_cons_of_mem _ ht))
        (fun t ht => hM t (List.mem_cons_of_mem _ ht))
        (fun t ht => hrow t (List.mem_cons_of_mem _ ht)) hf
      set C := B * ratioBound R rest with hC
      have hC0 : 0 ≤ C := by
        have := hrest x
        exact le_trans (abs_nonneg _) this
      have hinner : |∑ y, s.M x y * annealTest rest f y| ≤ C := by
        calc |∑ y, s.M x y * annealTest rest f y|
            ≤ ∑ y, |s.M x y * annealTest rest f y| := Finset.abs_sum_le_sum_abs _ _
          _ = ∑ y, s.M x y * |annealTest rest f y| := by
              refine Finset.sum_congr rfl (fun y _ => ?_)
              rw [abs_mul, abs_of_nonneg (hM s hs x y)]
          _ ≤ ∑ y, s.M x y * C := by
              refine Finset.sum_le_sum (fun y _ => ?_)
              exact mul_le_mul_of_nonneg_left (hrest y) (hM s hs x y)
          _ = C := by rw [← Finset.sum_mul, hrow s hs x, one_mul]
      show |(s.gnext x / s.gprev x) * ∑ y, s.M x y * annealTest rest f y| ≤ B * (R s * ratioBound R rest)
      obtain ⟨hr0, hr1⟩ := hr s hs x
      rw [abs_mul, abs_of_nonneg hr0]
      calc (s.gnext x / s.gprev x) * |∑ y, s.M x y * annealTest rest f y|
          ≤ R s * C := mul_le_mul hr1 hinner (abs_nonneg _) (hR s hs)
        _ = B * (R s * ratioBound R rest) := by rw [hC]; ring

/-- Bias of an inexact start. If the start is `μ₀` rather than `π₀`, the expected
weighted test function differs by at most `B * ∏ R_k * Σ_x |μ₀ x - π₀ x|`. -/
theorem bias_bound (R : AnnealStep Ω → ℝ) (steps : List (AnnealStep Ω)) (f : Ω → ℝ) (B : ℝ)
    (hR : ∀ s ∈ steps, 0 ≤ R s)
    (hr : ∀ s ∈ steps, ∀ x, 0 ≤ s.gnext x / s.gprev x ∧ s.gnext x / s.gprev x ≤ R s)
    (hM : ∀ s ∈ steps, ∀ x y, 0 ≤ s.M x y) (hrow : ∀ s ∈ steps, ∀ x, ∑ y, s.M x y = 1)
    (hf : ∀ x, |f x| ≤ B) (μ₀ π₀ : Ω → ℝ) :
    |∑ y, annealMeasure steps μ₀ y * f y - ∑ y, annealMeasure steps π₀ y * f y|
      ≤ B * ratioBound R steps * ∑ x, |μ₀ x - π₀ x| := by
  rw [annealMeasure_dual, annealMeasure_dual, ← Finset.sum_sub_distrib]
  have hb := annealTest_bound R steps f B hR hr hM hrow hf
  calc |∑ x, (μ₀ x * annealTest steps f x - π₀ x * annealTest steps f x)|
      ≤ ∑ x, |μ₀ x * annealTest steps f x - π₀ x * annealTest steps f x| := Finset.abs_sum_le_sum_abs _ _
    _ = ∑ x, |μ₀ x - π₀ x| * |annealTest steps f x| := by
        refine Finset.sum_congr rfl (fun x _ => ?_)
        rw [← sub_mul, abs_mul]
    _ ≤ ∑ x, |μ₀ x - π₀ x| * (B * ratioBound R steps) := by
        refine Finset.sum_le_sum (fun x _ => ?_)
        exact mul_le_mul_of_nonneg_left (hb x) (abs_nonneg _)
    _ = B * ratioBound R steps * ∑ x, |μ₀ x - π₀ x| := by
        rw [← Finset.sum_mul]; ring

end Tempering
