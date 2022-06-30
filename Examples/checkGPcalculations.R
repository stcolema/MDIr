
library(mdiHelpR)
library(tagmReDraft)
library(pheatmap)
library(ggplot2)
library(magrittr)

set.seed(1)

N <- 125
K <- 4
P <- 10
delta_mu <- 1.4
X <- generateSimulationDataset(K, N, P, delta_mu = delta_mu)
X$data <- X$data[ , findOrder(t(X$data))]
annotatedHeatmap(X$data, X$cluster_IDs, cluster_cols = F)

R <- 5000
thin <- 50
type <- "GP"
mcmc_out <- callMixtureModel(X$data, R, thin, type, K = 15)

burn <- 0.2 * R
discarded_samples <- seq(1, floor(burn / thin) + 1)

alloc <- mcmc_out$allocations[-discarded_samples, ]

psm <- makePSM(alloc)
row.names(psm) <- colnames(psm) <- row.names(X$data)
pheatmap(psm)
sim_col <- simColPal()
my_breaks <- defineBreaks(sim_col, lb = 0)
annotatedHeatmap(psm, X$cluster_IDs, my_breaks = my_breaks, col_pal = sim_col)
