#!/usr/bin/env Rscript

#' @title Generate simulation dataset
#' @description Generates a dataset based upon a mixture of $K$ Gaussian 
#' distributions with $p$ independent, relevant features and $p_n$ irrelevant 
#' features.
#' @param K The number of components to sample from.
#' @param n The number of samples to draw.
#' @param p The number of relevant (i.e. signal-bearing) features.
#' @param delta_mu The difference between the means defining each component 
#' within each feature (defaults to 1).
#' @param cluster_sd The standerd deviation of the Gaussian distributions.
#' @param pi The K-vector of the populations proportions across each component.
#' @param p_n The number of irrelevant features (defaults to 0).
#' @return A list of `data` (a data.frame of the generated data) and
#' `cluster_IDs` (a vector of the cluster membership of each item).
#' @export
generateSimulationDataset <- function(K, n, p,
                                delta_mu = 1,
                                cluster_sd = 1,
                                pi = rep(1/K, K),
                                p_n = 0) {
  
  # Create an empty list to hold the output
  my_data <- list(
    data = NA,
    cluster_IDs = NA
  )
  
  # Generate some cluster means and centre upon 0.
  cluster_means <- seq(from = 0, to = (K - 1) * delta_mu, by = delta_mu) %>%
    scale(center = T, scale = F)
  
  # If components overlap, crreate a K-vector of 0's
  if (delta_mu == 0) {
    cluster_means <- rep(0, K)
  }
  
  # Generate some cluster standard deviations
  # cluster_sds <- rep(cluster_sd, K)
  
  # if (pi_method == "even") {
  #   pi <- rep(1, K)
  # } else {
  #   pi <- rgamma(K, alpha)
  #   pi <- pi / sum(pi)
  # }
  
  # Find the number of requested informative features
  # p <- p # max(0, (p - p_noisy))
  
  # Define the global dataset standard deviation (this will be used to generate 
  # irrelevant features and will be updated after relevant features are 
  # generated)
  data_sd <- 1
  
  # Generate signal-bearing data if any relevant features are present
  if (p > 0) {
    my_data <- generateGaussianDataset(cluster_means, cluster_sd, n, p, pi)
    data_sd <- sd(my_data$data)
  }
  
  # If irrelevant features are desired, generate such data
  if (p_n > 0) {
    noisy_data <- lapply(1:p_n, function(x) {
      rnorm(n, sd = data_sd)
    }) %>%
      unlist() %>%
      matrix(ncol = p_n) %>%
      set_colnames(paste0("Noise_", 1:p_n))
    
    
    if (p > 0) {
      
      # Merge the relevant and irrelevant components of the data
      my_data$data <- cbind(my_data$data, noisy_data)
      
    } else {
      
      # If there are no relevant features make sure we handle this correctly
      my_data$data <- noisy_data %>%
        set_rownames(paste0("Person_", 1:n))
      
      # If no relevant features all items belong to the same population
      my_data$cluster_IDs <- rep(1, n)
      
    }
  }
  
  my_data
  
}
