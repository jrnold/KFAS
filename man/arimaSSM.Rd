\name{arimaSSM}
\alias{arimaSSM}
\title{Create a State Space Representation of ARIMA Model}
\usage{
  arimaSSM(y, arima, H = NULL, Q = NULL, u = NULL,
    distribution = c("Gaussian", "Poisson", "Binomial"),
    transform = c("none", "ldl", "augment"),
    tolF = .Machine$double.eps^0.5,
    tol0 = .Machine$double.eps^0.5)
}
\arguments{
  \item{arima}{A list or a list of lists with components
  \code{ar}, \code{ma} and \code{d}, giving the
  autoregression and moving average coefficients, and the
  degree of differencing for each series.  If arima is a
  single list, it is assumed that all \eqn{p} series have
  same ARIMA structure. Otherwise first sublist gives the
  ARIMA structure of the first series etc.}

  \item{H}{A \eqn{p \times p}{p*p} covariance matrix (or
  \eqn{p \times p \times n}{p*p*n} array in of time-varying
  case) of the disturbance terms
  \eqn{\epsilon_t}{\epsilon[t]} of the observation
  equation. Default gives \eqn{p \times p}{p*p} zero matrix
  ie. ordinary ARIMA model without additional noise.
  Omitted in case of non-Gaussian distributions. Augment
  the state vector if you to add want additional noise.}

  \item{Q}{A \eqn{p \times p}{p*p} covariance matrix of the
  disturbance terms \eqn{\eta_t}{\eta[t]} of the system
  equation. Default is \eqn{p \times p}{p*p} identity
  matrix ie. ordinary ARIMA model with disturbance terms
  having unit variance.}

  \item{y}{A time series object of class \code{ts}, or a
  object that can be coerced to such.}

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
\description{
  Function \code{arimaSSM} creates a state space
  representation of ARIMA model.
}
\details{
  The linear Gaussian state space model is given by

  \deqn{y_t = Z_t \alpha_t + \epsilon_t,}{y[t] =
  Z[t]\alpha[t] + \epsilon[t], (observation equation)}

  \deqn{\alpha_{t+1} = T_t \alpha_t + R_t
  \eta_t,}{\alpha[t+1] = T[t]\alpha[t] + R[t]\eta[t],
  (transition equation)}

  where \eqn{\epsilon_t ~ N(0,H_t)}{\epsilon[t] ~
  N(0,H[t])}, \eqn{\eta_t ~ N(0,Q_t)}{\eta[t] ~ N(0,Q[t])}
  and \eqn{\alpha_1 ~ N(a_1,P_1)}{\alpha[1] ~ N(a[1],P[1])}
  independently of each other. In case of non-Gaussian
  observations, the observation equation is of form
  \eqn{p(y_t|\theta_t) =
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
  log(\pi[t]/(1-\pi[t]))}, where \eqn{\pi_t}{\pi[t]} is the
  probability of success at time \eqn{t}.

  For non-Gaussian models \eqn{u_t=1}{u[t]=1} as a default.
  For Gaussian models, parameter is omitted.

  Only univariate observations are supported when
  observation equation is non-Gaussian.
}
\seealso{
  \code{\link{regSSM}} for state space representation of a
  regression model, \code{\link{structSSM}} for structural
  time series model, and \code{\link{SSModel}} for custom
  \code{SSModel} object.
}

