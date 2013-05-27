\name{ldl}
\alias{ldl}
\title{LDL Decomposition of a Matrix}
\usage{
  ldl(x, tol = .Machine$double.eps^0.5)
}
\arguments{
  \item{x}{Symmetrix matrix.}

  \item{tol}{Tolerance parameter for LDL decomposition,
  determines which diagonal values are counted as zero.}
}
\description{
  Function \code{ldl} computes the LDL decomposition of a
  symmetric matrix.
}

