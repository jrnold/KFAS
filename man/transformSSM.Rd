\name{transformSSM}
\alias{transformSSM}
\title{Transform the SSModel object with multivariate observations}
\usage{
  transformSSM(object, type = c("ldl", "augment"))
}
\arguments{
  \item{object}{State space model object from function
  SSModel.}

  \item{type}{Option \code{"ldl"} performs LDL
  decomposition for covariance matrix \eqn{H_t}{H[t]}, and
  multiplies the observation equation with the
  \eqn{L_t^{-1}}{L[t]^-1}, so \eqn{\epsilon_t^* \sim
  N(0,D_t)}{\epsilon[t]* ~ N(0,D[t])}. Option
  \code{"augment"} adds \eqn{\epsilon_t}{\epsilon[t]} to
  the state vector, when \eqn{Q_t}{Q[t]} becomes block
  diagonal with blocks \eqn{Q_t}{Q[t]} and \eqn{H_t}{H[t]}.
  In case of univariate series, option \code{"ldl"} only
  changes the \code{H_type} argument of the model to
  \code{"Diagonal"}.}
}
\value{
  \item{model}{Transformed model.}
}
\description{
  Function transform.SSModel transforms original model by
  LDL decomposition or state vector augmentation,
}
\details{
  As all the functions in KFAS use univariate approach,
  \eqn{H_t}{H[t]}, a covariance matrix of an observation
  equation needs to be either diagonal or zero matrix.
  Function transformSSM performs either the LDL
  decomposition of the covariance matrix of the observation
  equation, or augments the state vector with the
  disturbances of the observation equation.

  In case of a LDL decomposition, the new \eqn{H_t}{H[t]}
  contains the diagonal part of the decomposition, whereas
  observations \eqn{y_t}{Z[t]} and system matrices
  \eqn{Z_t}{Z[t]} are multiplied with the inverse of
  \eqn{L_t}{L[t]}.
}

