% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PUMICE.nested_cv.R
\name{elnet.cv}
\alias{elnet.cv}
\title{Nested Cross-validation}
\usage{
elnet.cv(
  geno_training,
  expression_training,
  fold_training,
  geno_validating,
  expression_validating,
  geno_tuning,
  expression_tuning,
  fold_tuning,
  geno_testing,
  expression_testing,
  penalty_k
)
}
\arguments{
\item{expression_training}{training expression set}

\item{fold_training}{index for traning for fold}

\item{geno_validating}{validating genotype set}

\item{expression_validating}{validating expression set}

\item{geno_tuning}{tuning genotype set}

\item{expression_tuning}{tuning expression set}

\item{fold_tuning}{index for tuning for fold}

\item{geno_testing}{testing genotype set}

\item{expression_testing}{testing expression set}

\item{penalty_k}{penalty factor}

\item{geno_traning}{training genotype set}
}
\description{
Nested Cross-validation
}
\examples{
elnet.cv()
}
\keyword{Nested}
\keyword{cross-validation}
