#!/usr/bin/env Rscript

#' @title Generate simulation dataset
#' @description Generates a dataset based upon a mixture of $K$ Gaussian
#' distributions with $P$ independent, relevant features and $P_n$ irrelevant
#' features. Irrelevant features contain no signal for underlying structure and
#' all measurements for an irrelevant feature are drawn from a common standard
#' Gaussian distribution.
#' @param K The number of components to sample from.
#' @param N The number of samples to draw.
#' @param P The number of relevant (i.e. signal-bearing) features.
#' @param delta_mu The difference between the means defining each component
#' within each feature (defaults to 1).
#' @param cluster_sd The standerd deviation of the Gaussian distributions.
#' @param pi The K-vector of the populations proportions across each component.
#' @param P_n The number of irrelevant features (defaults to 0).
#' @return A list of `data` (a data.frame of the generated data) and
#' `cluster_IDs` (a vector of the cluster membership of each item).
#' @importFrom stats sd
#' @returns Named list containing ``data``, a matrix of the generated Gaussian
#' data and ``cluster_IDs``, the true generating structure.
#' @examples
#' K <- 4
#' N <- 100
#' P <- 4
#' generateSimulationDataset(K, N, P)
#' @export
generateSimulationDataset <- function(K, N, P,
                                      delta_mu = 1,
                                      cluster_sd = 1,
                                      pi = rep(1 / K, K),
                                      P_n = 0) {

  # Create an empty list to hold the output
  my_data <- list(
    data = NA,
    cluster_IDs = NA
  )

  # Generate some cluster means and centre upon 0.
  cluster_means <- seq(from = 0, to = (K - 1) * delta_mu, by = delta_mu)
  cluster_means <- scale(cluster_means, center = TRUE, scale = FALSE)

  # If components overlap, crreate a K-vector of 0's
  if (delta_mu == 0) {
    cluster_means <- rep(0, K)
  }

  # Define the global dataset standard deviation (this will be used to generate
  # irrelevant features and will be updated after relevant features are
  # generated)
  data_sd <- 1

  # Generate signal-bearing data if any relevant features are present
  if (P > 0) {
    my_data <- generateGaussianDataset(cluster_means, cluster_sd, N, P, pi)
    data_sd <- stats::sd(my_data$data)
  }

  # If irrelevant features are desired, generate such data
  if (P_n > 0) {
    .noise <- lapply(seq(1, P_n), function(x) {
      rnorm(N, sd = data_sd)
    })
    noisy_data <- matrix(unlist(.noise), ncol = P_n)
    colnames(noisy_data) <- paste0("Noise_", seq(1, P_n))

    if (P > 0) {
      # Merge the relevant and irrelevant components of the data
      my_data$data <- cbind(my_data$data, noisy_data)
    } else {
      # If there are no relevant features make sure we handle this correctly
      my_data$data <- noisy_data
      row.names(my_data$data) <- paste0("Person_", seq(1, N))

      # If no relevant features all items belong to the same population
      my_data$cluster_IDs <- rep(1, N)
    }
  }

  my_data
}
