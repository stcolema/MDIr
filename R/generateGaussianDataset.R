#!/usr/bin/env Rscript

#' @title Generate Gaussian dataset
#' @description Generate a dataset based upon a mixture of Gaussian distributions
#' (with independent features).
#' @param cluster_means A k-vector of cluster means defining the k clusters.
#' @param n The number of samples to generate in the entire dataset.
#' @param p The number of columns to generate in the dataset.
#' @param pi A k-vector of the expected proportion of points to be drawn from
#' each distribution.
#' @param row_names The row names of the generated dataset.
#' @param col_names The column names of the generated dataset.
#' @export
generateGaussianDataset <- function(cluster_means, std_dev, n, p, pi,
                            row_names = paste0("Person_", 1:n),
                            col_names = paste0("Gene_", 1:p)) {
  
  # The number of distirbutions to sample from
  K <- length(cluster_means)
  
  # The membership vector for the n points
  cluster_IDs <- sample(K, n, replace = T, prob = pi)
  
  # The data matrix
  my_data <- matrix(nrow = n, ncol = p)
  
  # Sample the cluster mean order over the p columns
  # reordered_cluster_means <- lapply(1:p, function(x){sample(cluster_means)}) %>% 
  #   unlist() %>% 
  #   matrix(ncol = p, nrow = K, byrow = F)
  # 
  # # Create a N x P matrix of the cluster means associated with each item
  # component_struc <- reordered_cluster_means[cluster_IDs, ]
  # 
  # # Sample from the normal distributions associated with each item
  # my_data <- component_struc %>% 
  #   apply(2, function(x){MASS::mvrnorm(1, mu = x, diag(1, nrow = n))})
  # 
  # for(i in 1:n){
  #   my_data[i, ] <- reordered_cluster_means[cluster_IDs[i], ] %>%
  #     lapply(function(x){rnorm(p, mean = x, sd = std_dev)})

  # Iterate over each of the columns permuting the means associated with each
  # label.
  for (j in 1:p)
  {
    reordered_cluster_means <- sample(cluster_means)

    # Draw n points from the K univariate Gaussians defined by the permuted means.
    for (i in 1:n) {
      my_data[i, j] <- rnorm(1,
                             mean = reordered_cluster_means[cluster_IDs[i]],
                             sd = std_dev
      )
    }
  }
  
  # Order based upon allocation label
  row_order <- order(cluster_IDs)
  
  # Assign rownames and column names
  rownames(my_data) <- row_names
  colnames(my_data) <- col_names
  
  # Return the data and the allocation labels
  list(
    data = my_data[row_order, ],
    cluster_IDs = cluster_IDs[row_order]
  )
}