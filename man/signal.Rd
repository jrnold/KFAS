\name{signal}
\alias{signal}
\title{Extract the filtered or smoothed signal of a State Space Model}
\usage{
  signal(object, states = NULL)
}
\arguments{
  \item{object}{Object of class \code{KFS}.}

  \item{states}{Which states are combined? Default is NULL
  which combines all states according to \eqn{Z_t}{Z[t]}.}
}
\value{
  \item{signal}{Time series object of filtered signal
  \eqn{Z_ta_t}{Z[t]a[t]} or smoothed signal
  \eqn{Z_t\hat\alpha_t}{Z[t]\alpha[t]} using only the
  defined states. Filtered signal is computed only for
  non-diffuse phase.  }
  \item{variance}{Cov(\eqn{Z_ta_t}{Z[t]a[t]}) or
  Cov(\eqn{Z_t\hat\alpha_t}{Z[t]\alpha[t]}) using only the
  defined states. Covariances of filtered signal is
  computed only for non-diffuse phase.  }
}
\description{
  Function \code{signal} extracts the filtered or smoothed
  signal of a State Space model depending on the input
  object.
}

