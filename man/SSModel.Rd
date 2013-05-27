\name{SSModel}
\alias{SSModel}
\title{Create a State Space Model Object of class \code{SSModel}}
\usage{
  SSModel(y, Z = NULL, H = NULL, T = NULL, R = NULL,
    Q = NULL, a1 = NULL, P1 = NULL, P1inf = NULL, u = NULL,
    distribution = c("Gaussian", "Poisson", "Binomial"),
    transform = c("none", "ldl", "augment"),
    tolF = .Machine$double.eps^0.5,
    tol0 = .Machine$double.eps^0.5)
}
\arguments{
  \item{y}{A time series object of class \code{ts}, or a
  object that can be coerced to such.}

  \item{Z}{System matrix or array of observation equation.}

  \item{H}{Covariance matrix or array of disturbance terms
  \eqn{\epsilon_t}{\epsilon[t]} of observation equation.
  Omitted in case of non-Gaussian distributions. Augment
  the state vector if you want to add additional noise.}

  \item{T}{System matrix or array of transition equation.}

  \item{R}{System matrix or array of transition equation.}

  \item{Q}{Covariance matrix or array of disturbance terms
  \eqn{\eta_t}{\eta[t]}.}

  \item{a1}{Expected value of the initial state vector
  \eqn{\alpha_1}{\alpha[1]}.}

  \item{P1}{Covariance matrix of \eqn{\alpha_1}{\alpha[1]}.
  In the diffuse case the non-diffuse part of
  \eqn{P_1}{P[1]}.}

  \item{P1inf}{Diffuse part of \eqn{P_1}{P[1]}. Diagonal
  matrix with ones on diagonal elements which correspond to
  the unknown initial states.}

  \item{u}{Only used with non-Gaussian distribution. See
  details.}

  \item{distribution}{Specify the distribution of the
  observations. Default is "Gaussian".}

  \item{transform}{The functions of \code{KFAS} require
  diagonal covariance matrix \eqn{H_t}{H[t]}. If
  \eqn{H_t}{H[t]} is not diagonal, model can be transformed
  using one of the two options. Option \code{"ldl"}
  performs LDL decomposition for covariance matrix
  \eqn{H_t}{H[t]}, and multiplies the observation equation
  with the \eqn{L_t^{-1}}{L[t]^{-1}}, so \eqn{\epsilon_t
  \sim N(0,D_t)}{\epsilon[t] ~ N(0,D[t])}. Option
  \code{"augment"} adds \eqn{\epsilon_t}{\epsilon[t]} to
  the state vector, when \eqn{Q_t}{Q[t]} becomes block
  diagonal with blocks \eqn{Q_t}{Q[t]} and \eqn{H_t}{H[t]}.
  In case of univariate series, option \code{"ldl"} only
  changes the \code{H_type} argument of the model to
  \code{"Diagonal"}. Default is \code{"none"} which does no
  transformation but checks if \eqn{H} is diagonal. If not,
  \code{H_type} is set to \code{"Untransformed"}.}

  \item{tolF}{Tolerance parameter for Finf.  Smallest value
  not counted for zero.}

  \item{tol0}{Tolerance parameter for LDL decomposition,
  determines which diagonal values are counted as zero.}
}
\value{
  object of class \code{SSModel} with elements
}
\description{
  Function \code{SSModel} creates a state space object
  object of class \code{SSModel} which can be used as an
  input object for various functions of \code{KFAS}
  package.
}
\details{
  The custom state space model is constructed by using the
  given system matrices \code{Z}, \code{H}, \code{T},
  \code{R}, \code{Q}, \code{a1}, \code{P1} and
  \code{P1inf}. Matrix or scalar \code{Z} (array in case of
  time-varying \code{Z}) is used to determine the number of
  states \eqn{m}. If some of the other elements of the
  object are missing, \code{SSModel} uses default values
  which are identity matrix for \code{T}, \code{R} (or
  \eqn{k} first columns of identity matrix) and
  \code{P1inf}, and zero matrix for \code{H}, \code{Q},
  \code{P1} and , \code{a1}. If \code{P1} is given and
  \code{P1inf} is not, the it is assumed to be zero matrix.
  If \code{Q} is given, it is used to define \eqn{r}, the
  dimensions of \code{Q}, which can be smaller than \eqn{m}
  (defaults to \eqn{m}).

  The linear Gaussian state space model is given by

  \deqn{y_t = Z_t \alpha_t + \epsilon_t,}{y[t] =
  Z[t]\alpha[t] + \epsilon[t], (observation equation)}

  \deqn{\alpha_{t+1} = T_t \alpha_t + R_t
  \eta_t,}{\alpha[t+1] = T[t]\alpha[t] + R[t]\eta[t],
  (transition equation)}

  where \eqn{\epsilon_t \sim N(0,H_t)}{\epsilon[t] ~
  N(0,H[t])}, \eqn{\eta_t \sim N(0,Q_t)}{\eta[t] ~
  N(0,Q[t])} and \eqn{\alpha_1 \sim N(a_1,P_1)}{\alpha[1] ~
  N(a[1],P[1])} independently of each other. In case of
  non-Gaussian observations, the observation equation is of
  form \eqn{p(y_t|\theta_t) =
  p(y_t|Z_t\alpha_t)}{p(y[t]|\theta[t]) =
  p(y[t]|Z[t]\alpha[t])}, with
  \eqn{p(y_t|\theta_t)}{p(y[t]|\theta[t])} being one of the
  following:

  If observations are Poisson distributed, parameter of
  Poisson distribution is
  \eqn{u_t\lambda_t}{u[t]\lambda[t]} and \eqn{\theta_t =
  log(\lambda_t)}{\theta[t]=log(\lambda[t])}.

  If observations are from binomial distribution, \eqn{u}
  is a vector specifying number the of trials at times
  \eqn{1,\ldots,n}, and \eqn{\theta_t =
  log[\pi_t/(1-\pi_t)]}{\theta[t] =
  log(\pi[t]/(1-\pi[t]))}, where \eqn{\pi_t}{pi[t]} is the
  probability of success at time \eqn{t}.

  For non-Gaussian models \eqn{u_t=1}{u[t]=1} as a default.
  For Gaussian models, parameter is omitted.

  Only univariate observations are supported when
  observation equation is non-Gaussian.
}
\seealso{
  \code{\link{arimaSSM}} for state space representation of
  ARIMA model, \code{\link{regSSM}} for state space
  representation of a regression model, and
  \code{\link{structSSM}} for structural time series model.
}

