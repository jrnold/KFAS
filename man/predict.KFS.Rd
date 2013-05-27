\name{predict.KFS}
\alias{predict}
\alias{predict.KFS}
\title{Estimate the Missing Observations of a State Space Model}
\usage{
  \method{predict}{KFS} (object, fill = FALSE, ...)
}
\arguments{
  \item{object}{object of class \code{KFS} ie. the output
  from function KFS.}

  \item{fill}{If FALSE, only predictions of missing
  observations are returned, and other time points are
  markes as NA. This is convinient for plotting purposes.
  If TRUE, original time series is filled with predicted
  values for missing observations, and \eqn{F_t = 0}{F[t] =
  0} if observation is not missing. Default is FALSE.}

  \item{\dots}{Ignored.}
}
\value{
  A list with the following components: \item{y}{Time
  series object with missing observations replaced by
  \eqn{E(y_t|y)}{E(y[t]|y)}.  } \item{F}{Covariances
  \eqn{Cov(y_t|y)}{Cov(y[t]|y)}. Note that this is the
  usual multivariate version of \eqn{F_t}{F[t]} given by
  \eqn{Z_tP_tZ'_t + H_t}{Z[t]P[t]Z'[t] + H[t]}, not the
  univariate version given by KFS.  }
}
\description{
  Function \code{predict.KFS} computes the expected values
  of missing observations given the observed data.
}

