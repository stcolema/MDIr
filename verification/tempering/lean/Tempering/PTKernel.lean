import Mathlib
import Tempering.SwapKernel
import Tempering.ProductSwap

/-!
# One replica-exchange iteration is invariant

An iteration is: a local kernel at every temperature (independently), then a
round of exchanges between listed pairs of temperatures. For the deterministic
even-odd scheme the list is the even pairs on even rounds and the odd pairs on
odd rounds; for the stochastic scheme one of the two lists is chosen at random.
All we need is that every list of exchanges, and every choice of list, keeps
the product target invariant, which is proved here for arbitrary lists.
-/

open Finset

set_option linter.unusedSectionVars false
set_option linter.unusedSimpArgs false

namespace Tempering

variable {ι S : Type*} [Fintype ι] [DecidableEq ι] [Fintype S] [DecidableEq S]

/-- The kernel of a round of exchanges, applied in list order. -/
noncomputable def roundK (p : S → ℝ) (ℓ : S → ℝ) (β : ι → ℝ) :
    List (ι × ι) → (ι → S) → (ι → S) → ℝ
  | [] => fun s s' => if s' = s then 1 else 0
  | (ij :: rest) =>
      comp (swapK (prodTarget p ℓ β) (exch ij.1 ij.2)) (roundK p ℓ β rest)

/-- Any round of exchanges leaves the product target invariant, whatever the
list of pairs (so also for the even and odd rounds of DEO). -/
theorem roundK_invariant {p : S → ℝ} (hp : ∀ s, 0 < p s) (ℓ : S → ℝ) (β : ι → ℝ) :
    ∀ L : List (ι × ι), Invariant (prodTarget p ℓ β) (roundK p ℓ β L)
  | [] => invariant_id _
  | (ij :: rest) =>
      invariant_comp (exch_invariant hp ℓ β ij.1 ij.2) (roundK_invariant hp ℓ β rest)

/-- One iteration of parallel tempering: independent local kernels `P t`, each
invariant for its own tempered target, followed by a round of exchanges. The
product target `∏ π t` is invariant. -/
theorem pt_iteration_invariant {p : S → ℝ} (hp : ∀ s, 0 < p s) (ℓ : S → ℝ) (β : ι → ℝ)
    (P : ι → S → S → ℝ) (hP : ∀ t, Invariant (tw p ℓ β t) (P t)) (L : List (ι × ι)) :
    Invariant (prodTarget p ℓ β) (comp (prodKernel P) (roundK p ℓ β L)) := by
  apply invariant_comp
  · exact prodKernel_invariant (tw p ℓ β) P hP
  · exact roundK_invariant hp ℓ β L

end Tempering
