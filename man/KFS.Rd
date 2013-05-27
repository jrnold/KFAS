\name{KFS}
\alias{KFS}
\title{Kalman Filter and Smoother with Exact Diffuse Initialization for Exponential Family State Space Models}
\usage{
  KFS(object,
    smoothing = c("state", "disturbance", "both", "none"),
    simplify = TRUE, transform = c("ldl", "augment"),
    nsim = 100, theta = NULL, maxiter = 100)
}
\arguments{
  \item{object}{Object of class \code{SSModel} or
  \code{KFS} (in which case only smoothing is performed).}

  \item{smoothing}{Perform state or disturbance smoothing
  or both. Default is \code{"state"} for Gaussian models.
  For non-Gaussian models, state smoothing is always
  performed.}

  \item{simplify}{If FALSE, KFS returns some generally not
  so interesting variables from filtering and smoothing.
  Default is TRUE.}

  \item{transform}{How to transform the model in case of
  non-diagonal covariance matrix \eqn{H}. Defaults to
  \code{"ldl"}. See function \code{\link{transformSSM}} for
  details.}

  \item{nsim}{Number of independent samples. Default is
  100. Only used for non-Gaussian model.}

  \item{theta}{Initial values for conditional mode theta.
  Default is \code{log(mean(y/u))} for Poisson and
  \code{log(mean(y/(u-y)))} for Binomial distribution (or
  \code{log(mean(y))} in case of \eqn{u_t-y_t =
  0}{u[t]-y[t] = 0} for some \eqn{t}). Only used for
  non-Gaussian model.}

  \item{maxiter}{Maximum number of iterations used in
  linearisation. Default is 100. Only used for non-Gaussian
  model.}
}
\value{
  For Gaussian model, a list with the following components:
  \item{model}{Original state space model.  }
  \item{KFS.transform}{Type of H after possible
  transformation.  } \item{logLik}{Value of the
  log-likelihood function.  } \item{a}{One step predictions
  of states, \eqn{a_t=E(\alpha_t | y_{t-1}, \ldots ,
  y_{1})}{a[t]=E(\alpha[t] | y[t-1], \ldots , y[1])}.  }
  \item{P}{Covariance matrices of predicted states,
  \eqn{P_t=Cov(\alpha_t | y_{t-1}, \ldots ,
  y_{1})}{P[t]=Cov(\alpha[t] | y[t-1], \ldots , y[1])}.  }
  \item{Pinf}{Diffuse part of \eqn{P_t}{P[t]}. }
  \item{v}{Prediction errors \eqn{v_{i,t} = y_{i,t} -
  Z_{i,t}a_{i,t}, i=1,\ldots,p}{v[i,t] = y[i,t] -
  Z[i,t]a[i,t], i=1,\ldots,p},

  where \eqn{a_{i,t}=E(\alpha_t | y_{i-1,t}, \ldots,
  y_{1,t}, \ldots , y_{1,1})}{a[i,t]=E(\alpha[t] |
  y[i-1,t], \ldots, y[1,t], \ldots , y[1,1])}.  }

  \item{F}{Prediction error variances
  \eqn{Var(v_t)}{Var(v[t])}.  } \item{Finf}{Diffuse part of
  \eqn{F_t}{F[t]}.  } \item{d}{The last index of diffuse
  phase, i.e. the non-diffuse phase began from time
  \eqn{d+1}.  } \item{j}{The index of last \eqn{y_{i,t}} of
  diffuse phase.  } \item{alphahat}{Smoothed estimates of
  states, \eqn{E(\alpha_t | y_1, \ldots , y_n)}{E(\alpha[t]
  | y[1], \ldots , y[n])}. Only computed if
  \code{smoothing="state"} or \code{smoothing="both"}.  }
  \item{V}{Covariances \eqn{Var(\alpha_t | y_1, \ldots ,
  y_n).}{Var(\alpha[t] | y[1], \ldots , y[n]).} Only
  computed if \code{smoothing="state"} or
  \code{smoothing="both"}.  } \item{etahat}{Smoothed
  disturbance terms \eqn{E(\eta_t | y_1, \ldots ,
  y_n)}{E(\eta[t] | y[1], \ldots , y[n])}.Only computed if
  \code{smoothing="disturbance"} or
  \code{smoothing="both"}.  } \item{V_eta}{Covariances
  \eqn{Var(\eta_t | y_1, \ldots , y_n)}{Var(\eta[t] | y[1],
  \ldots , y[n])}. Only computed if
  \code{smoothing="disturbance"} or \code{smooth="both"}.
  } \item{epshat}{Smoothed disturbance terms
  \eqn{E(\epsilon_{t} | y_1, \ldots , y_n)}{E(\epsilon[t] |
  y[1], \ldots , y[n])}. Only computed if
  \code{smoothing="disturbance"} or
  \code{smoothing="both"}. } \item{V_eps}{Diagonal elements
  of \eqn{Var(\epsilon_{t} | y_1, \ldots ,
  y_n)}{Var(\epsilon[t] | y[1], \ldots , y[n])}. Note that
  due to the diagonalization, off-diagonal elements are
  zero. Only computed if \code{smoothing="disturbance"} or
  \code{smoothing="both"}.  } In addition, if argument
  \code{simplify=FALSE}, list contains following
  components: \item{K}{Covariances \eqn{Cov(\alpha_{t,i},
  y_{t,i} | y_{i-1,t}, \ldots, y_{1,t}, y_{t-1}, \ldots ,
  y_{1}), \quad i=1,\ldots,p}{Cov(\alpha[t,i], y[t,i] |
  y[i-1,t], \ldots, y[1,t], y[t-1], \ldots , y[1]),
  i=1,\ldots,p}.  } \item{Kinf}{Diffuse part of
  \eqn{K_t}{K[t]}.  } \item{r}{Weighted sums of innovations
  \eqn{v_{t+1}, \ldots , v_{n}}{v[t+1], \ldots , v[n]}.
  Notice that in literature t in \eqn{r_t}{r[t]} goes from
  \eqn{0, \ldots, n}. Here \eqn{t=1, \ldots, n+1}. Same
  applies to all r and N variables.  } \item{r0,
  r1}{Diffuse phase decomposition of \eqn{r_t}{r[t]}.  }
  \item{N}{Covariances \eqn{Var(r_t)}{Var(r[t])} .  }
  \item{N0, N1, N2}{Diffuse phase decomposition of
  \eqn{N_t}{N[t]}.  }

  For non-Gaussian model, a list with the following
  components: \item{model}{Original state space model with
  additional elements from function \code{approxSSM}.  }
  \item{alphahat}{Smoothed estimates of states
  \eqn{E(\alpha_t | y_1, \ldots , y_n)}{E(\alpha[t] | y[1],
  \ldots , y[n])}.  } \item{V}{Covariances
  \eqn{Var(\alpha_t | y_1, \ldots , y_n)}{Var(\alpha[t] |
  y[1], \ldots , y[n])}.  } \item{yhat}{A time series
  object containing smoothed means of observation
  distributions, with parameter
  \eqn{u_texp(\hat\theta_t)}{u[t]exp(thetahat[t])} for
  Poisson and
  \eqn{u_texp(\hat\theta_t)/(1+exp(\hat\theta_t))}{u[t]exp(thetahat[t])/(1+exp(thetahat[t])}.
  } \item{V.yhat}{a vector of length containing smoothed
  variances of observation distributions.  }
}
\description{
  Performs Kalman filtering and smoothing with exact
  diffuse initialization using univariate approach for
  exponential family state space models. For non-Gaussian
  models, state smoothing is provided with additional
  smoothed mean and variance of observations.
}
\details{
  Notice that in case of multivariate observations,
  \code{v}, \code{F}, \code{Finf}, \code{K} and \code{Kinf}
  are usually not the same as those calculated in usual
  multivariate Kalman filter. As filtering is done one
  observation element at the time, the elements of
  prediction error \eqn{v_t}{v[t]} are uncorrelated, and
  \code{F}, \code{Finf}, \code{K} and \code{Kinf} contain
  only the diagonal elemens of the corresponding covariance
  matrices.

  In rare cases of a very long diffuse initialization phase
  with highly correlated states, cumulative rounding errors
  in computing \code{Finf} and \code{Pinf} can sometimes
  cause the diffuse phase end too early. Changing the
  tolerance parameter \code{tolF} to smaller (or larger)
  should help.
}
\references{
  Koopman, S.J. and Durbin J. (2000).  Fast filtering and
  smoothing for non-stationary time series models, Journal
  of American Statistical Assosiation, 92, 1630-38.  \cr

  Koopman, S.J. and Durbin J. (2001).  Time Series Analysis
  by State Space Methods. Oxford: Oxford University Press.
  \cr

  Koopman, S.J. and Durbin J. (2003).  Filtering and
  smoothing of state vector for diffuse state space models,
  Journal of Time Series Analysis, Vol. 24, No. 1.  \cr
}

