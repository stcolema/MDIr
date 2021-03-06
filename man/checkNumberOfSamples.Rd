% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkNumberOfSamples.R
\name{checkNumberOfSamples}
\alias{checkNumberOfSamples}
\title{Check number of samples}
\usage{
checkNumberOfSamples(R, thin, verbose = FALSE)
}
\arguments{
\item{R}{The number of iterations in the sampler.}

\item{thin}{The factor by which the samples generated are thinned, e.g. if
``thin=50`` only every 50th sample is kept.}
}
\description{
Internal function to check the number of MCMC samples.
}
