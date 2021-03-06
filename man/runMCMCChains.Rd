% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runMCMCChains.R
\name{runMCMCChains}
\alias{runMCMCChains}
\title{Run MCMC Chains}
\usage{
runMCMCChains(
  X,
  n_chains,
  R,
  thin,
  initial_labels,
  fixed,
  types,
  K,
  alpha = NULL,
  initial_labels_as_intended = FALSE
)
}
\arguments{
\item{X}{Data to cluster. List of $L$ matrices with the $N$ items to cluster 
held in rows.}

\item{n_chains}{Integer. Number of MCMC chains to run.}

\item{R}{The number of iterations in the sampler.}

\item{thin}{The factor by which the samples generated are thinned, e.g. if
``thin=50`` only every 50th sample is kept.}

\item{initial_labels}{Initial clustering. $N x L$ matrix.}

\item{fixed}{Which items are fixed in their initial label. $N x L$ matrix.}

\item{alpha}{The concentration parameter for the stick-breaking prior and the
weights in the model.}

\item{initial_labels_as_intended}{Logical indicating if the passed initial 
labels are as intended or should ``generateInitialLabels`` be called.}

\item{type}{Character vector indicating density type to use. 'MVN'
(multivariate normal), 'TAGM' (t-adjust Gaussian mixture) or 'C' (categorical).}

\item{K_max}{Vector indicating the number of components to include (the upper
bound on the number of clusters in each dataset).}
}
\value{
A named list containing the sampled partitions, component weights and
phi parameters, model fit measures and some details on the model call.
}
\description{
Run multiple chains of Multiple Dataset Integration (MDI) using
the same inputs in each model run.
}
