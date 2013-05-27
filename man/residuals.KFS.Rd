\name{residuals.KFS}
\alias{residuals.KFS}
\title{Extract Standardized Residuals of Kalman Filter and Smoother output}
\usage{
  \method{residuals}{KFS} (object, ...)
}
\arguments{
  \item{object}{KFS object}

  \item{...}{Ignored.}
}
\description{
  Extract Standardized Residuals of Kalman Filter and
  Smoother output
}
\details{
  For object of class KFS, three types of residuals can be
  computed: Standardized one-step ahead prediction
  residuals are defined as \deqn{v_{i,t}/\sqrt{F_{i,t}},
  \quad i=1,\ldots,p,t=d+1,\ldots,n,}{v[i,t]F[i,t]^-0.5,
  i=1,\ldots,p, t=d+1,\ldots,n,} with residuals being
  undefined in diffuse phase.

  Residuals based on the smoothed disturbance terms are
  defined as \deqn{\hat
  \epsilon_{i,t}/\sqrt{Var(\epsilon_{i,t})}, \quad
  i=1,\ldots,p, t=1,\ldots,n,}{\hat
  \epsilon[i,t]Var(\epsilon[i,t])^-0.5, i=1,\ldots,p,
  t=1,\ldots,n,} and \deqn{L^{-1}_t \hat \eta_t, \quad
  t=1,\ldots,n,}{L^{-1}[t] \eta[t], t=1,\ldots,n,} where
  \eqn{L_t}{L[t]} is the lower triangular matrix from
  Cholesky decomposition of \eqn{V_{\eta,t}}{V[\eta,t]}
}

