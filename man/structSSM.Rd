\name{structSSM}
\alias{structSSM}
\title{Create a Structural Time Series State Space Model}
\usage{
  structSSM(y, trend = "level", seasonal = "none",
    X = NULL, H = NULL, Q.level = NULL, Q.slope = NULL,
    Q.seasonal = NULL, Q.regression = NULL, u = NULL,
    distribution = c("Gaussian", "Poisson", "Binomial"),
    transform = c("none", "ldl", "augment"),
    tolF = .Machine$double.eps^0.5,
    tol0 = .Machine$double.eps^0.5)
}
\arguments{
  \item{trend}{A character vector defining the type of the
  level component of the model. For multivariate series,
  either one type, it is assumed that all \eqn{p} series
  have same type of trend components. Possible values are
  \code{"level"} (local level model) and \code{"slope"}
  (local linear trend model). Default is \code{"level"}.}

  \item{seasonal}{A character vector defining the type of
  the seasonal component of the model. For multivariate
  series, it is assumed that all \eqn{p} series have same
  type of seasonal components.  Possible values are
  \code{"none"} (no seasonal), \code{"time"} (time domain
  form) and \code{"time"} (frequency domain form). The
  length of the seasonal pattern is taken as the frequency
  attribute of the time series object \code{y}. Default is
  \code{"none"}.}

  \item{X}{A \eqn{n \times k}{n*k} matrix of explanatory
  variables, with each column containing one explanatory
  variable. It is assumed that all \eqn{p} series use same
  explanatory variables.}

  \item{H}{A \eqn{p \times p}{p*p} covariance matrix (or
  \eqn{p \times p \times n}{p*p*n} array in time-varying
  case) of the disturbance terms
  \eqn{\epsilon_t}{\epsilon[t]} of the observation
  equation. Default gives \eqn{p \times p}{p*p} zero
  matrix.  Omitted in case of non-Gaussian distributions.
  Augment the state vector if you want to add additional
  noise.}

  \item{Q.level}{A scalar or \eqn{p \times p}{p*p}
  covariance matrix (or \eqn{p \times p \times n}{p*p*n}
  array in time-varying case) of the disturbance terms
  \eqn{\xi_t}{\xi[t]} corresponding to the level process
  \eqn{\mu_t}{\mu[t]}.  Default gives diagonal matrix with
  NA's on diagonal.}

  \item{Q.slope}{A scalar or \eqn{p \times p}{p*p}
  covariance matrix (or \eqn{p \times p \times n}{p*p*n}
  array in time-varying case) of the disturbance terms
  \eqn{\zeta_t}{\xi[t]} corresponding to the slope process
  \eqn{\nu_t}{\nu[t]}.  Default gives diagonal matrix with
  NA's on diagonal. Omitted if \code{trend="level"}.}

  \item{Q.seasonal}{scalar or A \eqn{p \times p}{p*p}
  covariance matrix (or \eqn{p \times p \times n}{p*p*n}
  array in time-varying case) of the disturbance terms
  \eqn{\omega_t}{\omega[t]} corresponding to the seasonal
  process \eqn{\gamma_t}{\gamma[t]}.  Default gives
  diagonal matrix with NA's on diagonal. Omitted if
  \code{seasonal="none"}.  There are several
  \eqn{\omega_t}{\omega[t]} processes in the frequency
  domain case, but they are all identically distributed, so
  only the (co)variance structure of one of them need to be
  defined.}

  \item{Q.regression}{A scalar or \eqn{xn \times xn}{xn*xn}
  covariance matrix (or \eqn{xn \times xn \times
  n}{xn*xn*n} array in time-varying case) of the
  disturbance terms corresponding to the regression
  coefficient processes.  Default gives zero matrix i.e.
  ordinary time-invariant regression.}

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
  Function \code{structSSM} creates a state space
  representation of structural time series.
}
\details{
  The structural time series model has form

  \deqn{y_t = \mu_t + \gamma_t + \epsilon_t, \quad
  \epsilon_t \sim N(0, H_t) }{y[t] = \mu[t] + \epsilon[t],
  \epsilon[t] ~ N(0, H[t])}

  \deqn{\mu_{t+1} = \mu_t + \nu_t + \xi_t, \quad \xi_t \sim
  N(0, Q_{level,t})}{\mu[t+1] = \mu[t] + \nu[t] + \xi[t],
  \xi[t] ~ N(0, Q[level,t])}

  \deqn{\nu_{t+1} = \nu_t + \zeta_t, \quad \zeta_t \sim
  N(0, Q_{slope,t})}{\nu[t+1] = \nu[t] + \zeta[t], \zeta[t]
  ~ N(0, Q[slope,t]) }

  with seasonal component being either time domain form
  \deqn{\gamma_{t+1} = -\sum_{j=1}^{s-1}\gamma_{t+1-j} +
  \omega_t, \quad \omega_t \sim N(0,Q_{seasonal,t}),
  }{\gamma[t+1] = -\gamma[t] - \ldots - \gamma[t-s+2] +
  \omega[t], \omega[t] ~ N(0,Q[seasonal,t]), }

  or frequency domain form where

  \deqn{\gamma_{t} = \sum_{j=1}^{\lfloor s/2
  \rfloor}\gamma_{j,t} }{\gamma[t] = \gamma[1,t] + \ldots +
  \gamma[[s/2],t], } \deqn{\gamma_{j,t+1} = \gamma_{j,t}
  cos\lambda_j + \gamma^{\ast}_{j,t} sin\lambda_j +
  \omega_{j,t}, }{\gamma_{j,t+1} = \gamma[j,t]
  cos\lambda[j] + \gamma*[j,t] sin\lambda[j] + \omega[j,t],
  } \deqn{\gamma^{\ast}_{j,t+1} = - \gamma_{j,t}
  sin\lambda_j + \gamma^{\ast}_{j,t} cos\lambda_j +
  \omega^\ast_{j,t}, j=1,\ldots, \lfloor s/2 \rfloor,
  }{\gamma*[j,t+1] = - \gamma[j,t] sin\lambda[j] +
  \gamma*[j,t] cos\lambda[j] + \omega*[j,t], j=1,...,
  [s/2], } with \eqn{\omega_{j,t}}{\omega[j,t]} and
  \eqn{\omega^\ast_{j,t}}{\omega*[j,t]} being independently
  distributed variables with \eqn{N(0,
  Q_{seasonal,t})}{N(0, Q[seasonal,t])} distribution and
  \eqn{\lambda_j = 2\pi j/s}{\lambda[j] = 2\pi j/s}.

  Explanatory variables can also be added to the model; in
  \code{structSSM} function it is assumed that same
  explanatory variables are used for all series. See
  \code{\link{regSSM}} and \code{\link{+}} for more
  complicated settings.
}
\seealso{
  \code{\link{arimaSSM}} for state space representation of
  ARIMA model, \code{\link{regSSM}} for state space
  representation of a regression model,
  \code{\link{SSModel}} for custom \code{SSModel} object
  and \code{\link{KFAS}} for general information regarding
  the package and examples of its usage.
}

