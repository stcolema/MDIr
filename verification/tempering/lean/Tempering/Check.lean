import Tempering.SwapKernel
import Tempering.ProductSwap
import Tempering.Conjugate
import Tempering.PTKernel

-- Every theorem should depend only on the standard axioms (no `sorryAx`).
#print axioms Tempering.swapK_detailedBalance
#print axioms Tempering.swapK_invariant
#print axioms Tempering.invariant_comp
#print axioms Tempering.invariant_mix
#print axioms Tempering.prodTarget_exch_ratio
#print axioms Tempering.acc_exch_eq
#print axioms Tempering.exch_invariant
#print axioms Tempering.prodKernel_invariant
#print axioms Tempering.marginal_of_product
#print axioms Tempering.gibbsK_invariant
#print axioms Tempering.tempered_square
#print axioms Tempering.tempered_power
#print axioms Tempering.roundK_invariant
#print axioms Tempering.pt_iteration_invariant
