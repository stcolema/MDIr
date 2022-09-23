#!/usr/bin/env Rscript
#' @title Generate Gaussian dataset
#' @description Generate a dataset based upon a mixture of Gaussian distributions
#' (with independent features).
#' @param cluster_means A k-vector of cluster means defining the k clusters.
#' @param std_dev A k-vector of cluster standard deviations defining the k clusters.
#' @param N The number of samples to generate in the entire dataset.
#' @param P The number of columns to generate in the dataset.
#' @param pi A k-vector of the expected proportion of points to be drawn from
#' each distribution.
#' @param row_names The row names of the generated dataset.
#' @param col_names The column names of the generated dataset.
#' @returns Named list of ``data``, the generated matrix and ``cluster_IDs``,
#' the generating structure.
#' @importFrom stats rnorm
#' @examples
#' cluster_means <- c(-2, 0, 2)
#' std_dev <- c(1, 1, 1.25)
#' N <- 100
#' P <- 5
#' pi <- c(0.3, 0.3, 0.4)
#' generateGaussianDataset(cluster_means, std_dev, N, P, pi)
#' @export
generateGaussianDataset <- function(cluster_means, std_dev, N, P, pi,
                                    row_names = paste0("Person_", seq(1, N)),
                                    col_names = paste0("Gene_", seq(1, P))) {

  # The number of distirbutions to sample from
  K <- length(cluster_means)

  # The membership vector for the N points
  cluster_IDs <- sample(K, N, replace = T, prob = pi)

  # The data matrix
  my_data <- matrix(nrow = N, ncol = P)

  # Iterate over each of the columns permuting the means associated with each
  # label.
  for (j in seq(1, P))
  {
    reordered_cluster_means <- sample(cluster_means)

    # Draw N points from the K univariate Gaussians defined by the permuted means.
    for (i in seq(1, N)) {
      my_data[i, j] <- stats::rnorm(1,
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
