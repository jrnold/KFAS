\name{importanceSSM}
\alias{importanceSSM}
\title{Importance Sampling of Non-Gaussian State Space Model}
\usage{
  importanceSSM(model, nsim = 1000, save.model = FALSE,
    theta = NULL, antithetics = TRUE, maxiter = 100)
}
\arguments{
  \item{model}{Non-Gaussian state space model object of
  class \code{SSModel}.}

  \item{nsim}{Number of independent samples. Default is
  1000.}

  \item{save.model}{Return the original model with the
  samples. Default is FALSE.}

  \item{theta}{Initial values for conditional mode theta.
  Default is \code{log(mean(y/u))} for Poisson and
  \code{log(mean(y/(u-y)))} for Binomial distribution (or
  \code{log(mean(y))} in case of \eqn{u_t-y_t =
  0}{u[t]-y[t] = 0} for some \eqn{t}).}

  \item{antithetics}{Logical. If TRUE, two antithetic
  variables are used in simulations, one for location and
  another for scale. Default is TRUE.}

  \item{maxiter}{Maximum number of iterations used in
  linearisation. Default is 100.}
}
\description{
  Importance Sampling of Non-Gaussian State Space Model.
}
\details{
  Function \code{importanceSSM} simulates states of the
  non-Gaussian state space model conditioned with the
  observations, returning the simulated samples of the
  states with the importance weights.
}

