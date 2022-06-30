

library(tagmReDraft)
library(mdiHelpR)
library(magrittr)

set.seed(1)

N <- 100
P <- 4
K <- 4
delta_mu <- 2

my_data <- generateSimulationDataset(K, N, P, delta_mu = delta_mu)

annotatedHeatmap(my_data$data, my_data$cluster_IDs)

R <- 1000
thin <- 25
type <- "MVN"
K_max <- 12
samples <- callMixtureModel(my_data$data, R, thin, type, K_max)