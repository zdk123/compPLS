% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PLSDA.R, R/SPLSDA.R
\name{plsDA_main}
\alias{plsDA_main}
\alias{splsDA_main}
\title{plsDA_main partial least squares discriminant analysis}
\usage{
plsDA_main(
  x,
  grouping,
  K,
  usePriors = FALSE,
  fold = 5,
  nboots = 999,
  n.core = 4,
  noise = 0,
  ...
)

splsDA_main(
  x,
  grouping,
  eta,
  K,
  usePriors = FALSE,
  fold = 5,
  nboots = 999,
  n.core = 4,
  noise = 0,
  ...
)
}
\arguments{
\item{x}{data with samples in rows, features are columns (not necessarily compositional x)}

\item{grouping}{a numeric vector or factor with sample classes (length should equal \code{nrow(x)})}

\item{K}{numeric vector containing number of components in the PLS model}

\item{usePriors}{use priors for very biased sample size between groups (ie - put strong penalty on misclassifying small groups)}

\item{fold}{number of partitions to randomly subsample for cross-validation}

\item{nboots}{number of bootstraps/permutations for estimating coefficient p-vals}

\item{n.core}{number of cores for paralellization of bootstraps}

\item{noise}{for very sparse components, some subsamples may have zero variance. Optionally, add some Gaussian noise to to avoid PLS errors}

\item{...}{additional arguments passed to plsDA}
}
\value{
a \code{plsDA} object that contains: the plsda model/object, \code{pvals}, the original data, \code{x}, and \code{groupings}

a \code{plsDA} object that contains: the plsda model/object, \code{pvals}, the original data, \code{x}, and \code{groupings}
}
\description{
The main wrapper for full Partial Least Squares discriminant analysis,
performing cross-validation to tune model parameters (here, number of components)
and do permutation tests (ie bootstrapping) to get pseudo-pvals estimates for model coefficients

The main wrapper for full sparse Partial Least Squares discriminant analysis,
performing cross-validation to tune model parameters (here, number of components)
and do permutation tests (ie bootstrapping) to get pseudo-pvals estimates for model coefficients
}
\seealso{
\code{\link{plsDA}}

\code{\link{plsDA}}
}
