\name{plot.KFS}
\alias{plot.KFS}
\title{Plot Ouput of Kalman Filter and Smoother}
\usage{
  \method{plot}{KFS} (x, y = NULL, type = "signal", ...)
}
\arguments{
  \item{x}{object of class \code{KFS}}

  \item{y}{Ignored.}

  \item{type}{Draw \code{"signals"} together with
  observations, or figures of \code{"state"} components.
  Default is \code{"signal"}. If smoothed values of states
  are available, those are used. Otherwise filtered values
  are used.}

  \item{...}{Graphical parameters for function plot.ts.}
}
\description{
  Plot Ouput of Kalman Filter and Smoother
}

