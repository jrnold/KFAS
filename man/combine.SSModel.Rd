\name{+.SSModel}
\alias{+.SSModel}
\title{Combine State Space Model Objects of class \code{SSModel}}
\usage{
  \method{+}{SSModel} (e1, e2)
}
\arguments{
  \item{e1}{,}

  \item{e2}{Models to be combined.}
}
\value{
  \item{model}{Combined model.}
}
\description{
  Second model needs to have either only duplicate time
  series with first model, or no identical series at all.
}

