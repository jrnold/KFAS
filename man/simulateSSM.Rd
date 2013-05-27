\name{simulateSSM}
\alias{simulateSSM}
\title{Simulation of a Gaussian State Space Model}
\usage{
  simulateSSM(object,
    sim = c("states", "disturbances", "observations"),
    nsim = 1, antithetics = FALSE, conditional = TRUE)
}
\arguments{
  \item{object}{Gaussian state space object.}

  \item{sim}{What to simulate. Note that all the
  simulations are done independently.}

  \item{nsim}{Number of independent samples. Default is 1.}

  \item{antithetics}{Use antithetic variables in
  simulation. Default is FALSE.}

  \item{conditional}{Simulations are conditional to data.
  If FALSE, the initial state \eqn{\alpha_1}{\alpha[1]} is
  set to \eqn{\hat \alpha_1}{alphahat[1]} computed by
  \code{KFS}, and all the observations are removed from the
  model. Default is TRUE.}
}
\description{
  Function \code{simulateSMM} simulates states,
  disturbances or missing observations of the Gaussian
  state space object conditionally on the data.
}
\details{
  Simulation smoother algorithm is from article by J.
  Durbin and S.J. Koopman (2002).

  Function can use two antithetic variables, one for
  location and other for scale, so output contains four
  blocks of simulated values which correlate which each
  other (ith block correlates negatively with (i+1)th
  block, and positively with (i+2)th block etc.).
}
\references{
  Durbin J. and Koopman, S.J. (2002). A simple and
  efficient simulation smoother for state space time series
  analysis, Biometrika, Volume 89, Issue 3
}

