import Mathlib
import Tempering.SwapKernel

/-!
# Tempered product target, replica exchange, and the cold marginal

State space `ι → S` (one state per temperature index, `ι` finite). With base
weight `p : S → ℝ` (prior, strictly positive), data log-likelihood `ℓ : S → ℝ`
and inverse temperatures `β : ι → ℝ`, replica `t` targets the unnormalised
density `π t s = exp (β t * ℓ s) * p s`, and the product target is
`Π s = ∏ t, π t (s t)`.

We prove, for `i ≠ j` and the exchange `σ s = s ∘ swap i j`:

* the ratio `Π (σ s) / Π s` is exactly `exp ((β i - β j) * (ℓ (s j) - ℓ (s i)))`
  (the prior `p` cancels), so the acceptance probability implemented in the
  sampler is the Metropolis acceptance of this proposal;
* hence the exchange kernel is `Π`-reversible and `Π`-invariant
  (from `Tempering.SwapKernel`);
* a product of per-temperature invariant kernels is `Π`-invariant;
* if each `π t` is normalised, the `t`-th marginal of `Π` is `π t`
  (in particular the cold replica has the posterior as its marginal).
-/

open Finset

set_option linter.unusedSectionVars false
set_option linter.unusedSimpArgs false

namespace Tempering

variable {ι S : Type*} [Fintype ι] [DecidableEq ι] [Fintype S] [DecidableEq S]

/-- A product over a finite index set splits off two distinct coordinates. -/
lemma prod_split (f : ι → ℝ) {i j : ι} (hij : i ≠ j) :
    ∏ t, f t = f i * f j * ∏ t ∈ (Finset.univ.erase i).erase j, f t := by
  rw [← Finset.mul_prod_erase Finset.univ f (Finset.mem_univ i)]
  have hj : j ∈ Finset.univ.erase i := Finset.mem_erase.mpr ⟨hij.symm, Finset.mem_univ j⟩
  rw [← Finset.mul_prod_erase (Finset.univ.erase i) f hj]
  ring

/-- The tempered weight of one replica. -/
noncomputable def tw (p : S → ℝ) (ℓ : S → ℝ) (β : ι → ℝ) (t : ι) (s : S) : ℝ :=
  Real.exp (β t * ℓ s) * p s

/-- The product target. -/
noncomputable def prodTarget (p : S → ℝ) (ℓ : S → ℝ) (β : ι → ℝ) (s : ι → S) : ℝ :=
  ∏ t, tw p ℓ β t (s t)

/-- Exchange of the states at temperatures `i` and `j`. -/
def exch (i j : ι) (s : ι → S) : ι → S := fun t => s (Equiv.swap i j t)

lemma exch_involutive (i j : ι) (s : ι → S) : exch i j (exch i j s) = s := by
  funext t; simp [exch]

lemma prodTarget_pos {p : S → ℝ} (hp : ∀ s, 0 < p s) (ℓ : S → ℝ) (β : ι → ℝ) (s : ι → S) :
    0 < prodTarget p ℓ β s :=
  Finset.prod_pos (fun _ _ => mul_pos (Real.exp_pos _) (hp _))

/-- Cross-multiplied form of the exchange ratio: the factors away from `i`, `j`
are unchanged by the exchange. -/
lemma prodTarget_exch_mul {p : S → ℝ} (ℓ : S → ℝ) (β : ι → ℝ) {i j : ι} (hij : i ≠ j)
    (s : ι → S) :
    prodTarget p ℓ β (exch i j s) * (tw p ℓ β i (s i) * tw p ℓ β j (s j))
      = prodTarget p ℓ β s * (tw p ℓ β i (s j) * tw p ℓ β j (s i)) := by
  unfold prodTarget
  rw [prod_split (fun t => tw p ℓ β t (exch i j s t)) hij, prod_split (fun t => tw p ℓ β t (s t)) hij]
  have hrest : ∏ t ∈ (Finset.univ.erase i).erase j, tw p ℓ β t (exch i j s t)
      = ∏ t ∈ (Finset.univ.erase i).erase j, tw p ℓ β t (s t) := by
    refine Finset.prod_congr rfl (fun t ht => ?_)
    have h1 : t ≠ j := (Finset.mem_erase.mp ht).1
    have h2 : t ≠ i := (Finset.mem_erase.mp (Finset.mem_erase.mp ht).2).1
    simp [exch, Equiv.swap_apply_of_ne_of_ne h2 h1]
  rw [hrest]
  have hi : exch i j s i = s j := by simp [exch]
  have hj : exch i j s j = s i := by simp [exch]
  simp only [hi, hj]
  ring

/-- The exchange ratio is the likelihood-only expression: the priors cancel. -/
theorem prodTarget_exch_ratio {p : S → ℝ} (hp : ∀ s, 0 < p s) (ℓ : S → ℝ) (β : ι → ℝ)
    {i j : ι} (hij : i ≠ j) (s : ι → S) :
    prodTarget p ℓ β (exch i j s) / prodTarget p ℓ β s
      = Real.exp ((β i - β j) * (ℓ (s j) - ℓ (s i))) := by
  have hpos := prodTarget_pos hp ℓ β s
  rw [div_eq_iff hpos.ne']
  have hw : 0 < tw p ℓ β i (s i) * tw p ℓ β j (s j) :=
    mul_pos (mul_pos (Real.exp_pos _) (hp _)) (mul_pos (Real.exp_pos _) (hp _))
  have key := prodTarget_exch_mul (p := p) ℓ β hij s
  have hratio : tw p ℓ β i (s j) * tw p ℓ β j (s i)
      = Real.exp ((β i - β j) * (ℓ (s j) - ℓ (s i))) * (tw p ℓ β i (s i) * tw p ℓ β j (s j)) := by
    unfold tw
    have : Real.exp (β i * ℓ (s j)) * Real.exp (β j * ℓ (s i))
        = Real.exp ((β i - β j) * (ℓ (s j) - ℓ (s i)))
          * (Real.exp (β i * ℓ (s i)) * Real.exp (β j * ℓ (s j))) := by
      rw [← Real.exp_add, ← Real.exp_add, ← Real.exp_add]
      congr 1; ring
    calc Real.exp (β i * ℓ (s j)) * p (s j) * (Real.exp (β j * ℓ (s i)) * p (s i))
        = (Real.exp (β i * ℓ (s j)) * Real.exp (β j * ℓ (s i))) * (p (s j) * p (s i)) := by ring
      _ = Real.exp ((β i - β j) * (ℓ (s j) - ℓ (s i)))
          * (Real.exp (β i * ℓ (s i)) * Real.exp (β j * ℓ (s j))) * (p (s j) * p (s i)) := by
        rw [this]
      _ = _ := by ring
  rw [hratio] at key
  have := key
  nlinarith [this, hw]

/-- The acceptance probability `min 1 (Π(σ s)/Π s)` equals
`exp (min 0 ((β i - β j) * (ℓ (s j) - ℓ (s i))))`, the form used by Syed et al. (2022,
eq. 6) with potential `V = -ℓ`. -/
theorem acc_exch_eq {p : S → ℝ} (hp : ∀ s, 0 < p s) (ℓ : S → ℝ) (β : ι → ℝ)
    {i j : ι} (hij : i ≠ j) (s : ι → S) :
    acc (prodTarget p ℓ β) (exch i j) s
      = Real.exp (min 0 ((β i - β j) * (ℓ (s j) - ℓ (s i)))) := by
  unfold acc
  rw [prodTarget_exch_ratio hp ℓ β hij s]
  rcases le_total 0 ((β i - β j) * (ℓ (s j) - ℓ (s i))) with h | h
  · rw [min_eq_left h, Real.exp_zero]
    exact min_eq_left (Real.one_le_exp h)
  · rw [min_eq_right h]
    exact min_eq_right (Real.exp_le_one_iff.mpr h)

/-- The exchange kernel of temperatures `i ≠ j` leaves the product target invariant. -/
theorem exch_invariant {p : S → ℝ} (hp : ∀ s, 0 < p s) (ℓ : S → ℝ) (β : ι → ℝ) (i j : ι) :
    Invariant (prodTarget p ℓ β) (swapK (prodTarget p ℓ β) (exch i j)) :=
  swapK_invariant _ (prodTarget_pos hp ℓ β) _ (exch_involutive i j)

/-- A product of per-temperature kernels. -/
noncomputable def prodKernel (P : ι → S → S → ℝ) (s s' : ι → S) : ℝ := ∏ t, P t (s t) (s' t)

/-- If `P t` leaves `π t` invariant for each `t`, the product kernel leaves the
product `∏ π t` invariant. -/
theorem prodKernel_invariant (π : ι → S → ℝ) (P : ι → S → S → ℝ)
    (h : ∀ t, Invariant (π t) (P t)) :
    Invariant (fun s : ι → S => ∏ t, π t (s t)) (prodKernel P) := by
  intro s'
  unfold prodKernel
  have e : ∀ s : ι → S, (∏ t, π t (s t)) * ∏ t, P t (s t) (s' t) = ∏ t, (π t (s t) * P t (s t) (s' t)) :=
    fun s => (Finset.prod_mul_distrib).symm
  simp_rw [e]
  have := Finset.prod_univ_sum (fun _ : ι => (Finset.univ : Finset S))
    (fun t (a : S) => π t a * P t a (s' t))
  rw [Fintype.piFinset_univ] at this
  rw [← this]
  exact Finset.prod_congr rfl (fun t _ => h t (s' t))

/-- If each `π t` sums to one, the `t₀` marginal of the product target is `π t₀`. -/
theorem marginal_of_product (π : ι → S → ℝ) (hnorm : ∀ t, ∑ a, π t a = 1) (t₀ : ι) (x : S) :
    ∑ s : ι → S, (if s t₀ = x then ∏ t, π t (s t) else 0) = π t₀ x := by
  have key := Finset.prod_univ_sum (fun _ : ι => (Finset.univ : Finset S))
    (fun t (a : S) => if t = t₀ then (if a = x then π t a else 0) else π t a)
  rw [Fintype.piFinset_univ] at key
  have lhs : ∀ s : ι → S, (if s t₀ = x then ∏ t, π t (s t) else 0)
      = ∏ t, (if t = t₀ then (if s t = x then π t (s t) else 0) else π t (s t)) := by
    intro s
    by_cases hs : s t₀ = x
    · rw [if_pos hs]
      refine Finset.prod_congr rfl (fun t _ => ?_)
      by_cases ht : t = t₀
      · subst ht; simp [hs]
      · simp [ht]
    · rw [if_neg hs]
      symm
      apply Finset.prod_eq_zero (Finset.mem_univ t₀)
      simp [hs]
  simp_rw [lhs]
  rw [← key]
  rw [← Finset.mul_prod_erase Finset.univ _ (Finset.mem_univ t₀)]
  have h1 : ∑ a, (if t₀ = t₀ then (if a = x then π t₀ a else 0) else π t₀ a) = π t₀ x := by
    simp
  have h2 : ∏ t ∈ Finset.univ.erase t₀, ∑ a, (if t = t₀ then (if a = x then π t a else 0) else π t a) = 1 := by
    refine Finset.prod_eq_one (fun t ht => ?_)
    have : t ≠ t₀ := (Finset.mem_erase.mp ht).1
    simp [this, hnorm t]
  rw [h1, h2, mul_one]

end Tempering
