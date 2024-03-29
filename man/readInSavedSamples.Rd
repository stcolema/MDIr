% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readInSavedSamples.R
\name{readInSavedSamples}
\alias{readInSavedSamples}
\title{Read in Saved Samples}
\usage{
readInSavedSamples(run_details)
}
\arguments{
\item{run_details}{Output from `callMDIWritingToFile` - a list describing the
model run, saved sample location, etc.}
}
\value{
A named list containing the sampled partitions, component weights,
phi and mass parameters and some details on the model call.
}
\description{
Reads in the saved output of `callMDIWritingToFile` and compiles
it into a more easily usable format.
}
\examples{

N <- 100
X <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)
Y <- matrix(c(rnorm(N, 0, 1), rnorm(N, 3, 1)), ncol = 2, byrow = TRUE)

truth <- c(rep(1, N / 2), rep(2, N / 2))
data_modelled <- list(X, Y)

V <- length(data_modelled)

# This R is much too low for real applications
R <- 100
thin <- 5
burn <- 10

K_max <- 10
K <- rep(K_max, V)
types <- rep("G", V)

model_call <- callMDIWritingToFile(data_modelled, R, thin, types, K = K)
mcmc_out <- readInSavedSamples(model_call)
}
