\name{logLik.SSModel}
\alias{logLik}
\alias{logLik.SSModel}
\title{Log-likelihood of the State Space Model.}
\usage{
  \method{logLik}{SSModel} (object, nsim = 0,
    antithetics = TRUE, taylor = TRUE, theta = NULL,
    maxiter = 100, fix.seed = TRUE, ...)
}
\arguments{
  \item{object}{State space model of class \code{SSModel}.}

  \item{nsim}{Number of independent samples used in
  estimating the log-likelihood of the non-gaussian state
  space model. Default is 0, which gives good starting
  value for optimization. Only used for non-Gaussian
  model.}

  \item{antithetics}{Logical. If TRUE, two antithetic
  variables are used in simulations, one for location and
  another for scale. Default is TRUE. Only used for
  non-Gaussian model.}

  \item{taylor}{Logical. If TRUE, control variable based on
  Taylor series is used. Default is TRUE. Only used for
  non-Gaussian model.}

  \item{theta}{Initial values for conditional mode theta.
  Default is \code{log(mean(y/u))} for Poisson and
  \code{log(mean(y/(u-y)))} for Binomial distribution (or
  \code{log(mean(y))} in case of \eqn{u_t-y_t =
  0}{u[t]-y[t] = 0} for some \eqn{t}). Only used for
  non-Gaussian model.}

  \item{fix.seed}{Use fixed seed. If FALSE, no fixed seed
  is used. If fix.seed is positive value, the value is used
  as a seed via set.seed function. Default is TRUE, so that
  the variation in random number generation does not affect
  numerical optimization algorithms. Only used for
  non-Gaussian model.}

  \item{...}{Ignored.}

  \item{maxiter}{Maximum number of iterations used in
  linearisation. Default is 100.}
}
\value{
  \item{}{log-likelihood of the state space model.}
}
\description{
  Function \code{logLik.SSmodel} computes the
  log-likelihood value of a state-space model.
}

