\name{approxSSM}
\alias{approxSSM}
\title{Linear Gaussian Approximation for Non-Gaussian State Space Model}
\usage{
  approxSSM(object, theta = NULL, maxiter = 100)
}
\arguments{
  \item{object}{Non-Gaussian state space model object of
  class \code{SSModel}.}

  \item{theta}{Initial values for conditional mode theta.
  Default is \code{log(mean(y/u))} for Poisson and
  \code{log(mean(y/(u-y)))} for Binomial distribution (or
  \code{log(mean(y))} in case of \eqn{u_t-y_t =
  0}{u[t]-y[t] = 0} for some \eqn{t}).}

  \item{maxiter}{Maximum number of iterations used in
  linearisation. Default is 100.}
}
\value{
  An object which contains the approximating Gaussian state
  space model with additional components
  \code{original.distribution}, \code{original.y},
  \code{thetahat}, and \code{iterations} (the number of
  iterations used).
}
\description{
  Function \code{approxSMM} computes the linear Gaussian
  approximation of a state space model where the
  observations have a non-Gaussian exponential family
  distribution. Currently only Poisson and Binomial
  distributions are supported.
}
\details{
  The linear Gaussian approximating model is a model
  defined by \deqn{\tilde y_t = Z_t \alpha_t + \epsilon_t,
  \quad \epsilon_t \sim N(0,\tilde H_t),}{ytilde[t] =
  Z[t]\alpha[t] + \epsilon[t], \epsilon[t] ~
  N(0,Htilde[t]),} \deqn{\alpha_{t+1} = T_t \alpha_t + R_t
  \eta_t, \quad \eta_t \sim N(0,Q_t),}{\alpha[t+1] =
  T[t]\alpha[t] + R[t]\eta[t], \eta[t] ~ N(0,Q[t]),} and
  \eqn{\alpha_1 \sim N(a_1,P_1)}{\alpha[1] ~ N(a[1],P[1])},
  where \eqn{\tilde y}{ytilde} and \eqn{\tilde H}{Htilde}
  is chosen in a way that the linear Gaussian approximating
  model has the same conditional mode of
  \eqn{\theta=Z\alpha} given the observations \eqn{y} as
  the original non-Gaussian model. Models also have same
  curvature at the mode.

  The linearization of the exponential family state space
  model is based on the first two derivatives of the
  observational logdensity.

  The approximating Gaussian model is used in computation
  of the log-likelihood of the non-Gaussian model and in
  importance sampling of non-Gaussian model.
}
\seealso{
  Importance sampling of non-Gaussian state space models
  \code{\link{importanceSSM}}, construct a \code{SSModel}
  object \code{\link{SSModel}}.
}

