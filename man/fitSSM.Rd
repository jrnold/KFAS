\name{fitSSM}
\alias{fitSSM}
\title{Maximum Likelihood Estimation of a State Space Model}
\usage{
  fitSSM(inits, model = NULL, modFun = NULL,
    method = "BFGS", nsim = 0, antithetics = TRUE,
    taylor = TRUE, theta = NULL, maxiter = 500, ...)
}
\arguments{
  \item{inits}{Initial values for \code{optim}}

  \item{model}{Model object of class \code{SSModel}. if
  \code{ModFun} is defined, this argument is ignored.}

  \item{modFun}{User defined function which builds the
  model of class \code{SSModel} given the parameters. If
  NULL, default estimation procedure is used (See
  details).}

  \item{method}{The method to be used in \code{optim}.
  Default is \code{"BFGS"}.}

  \item{nsim}{Number of independent samples used in
  estimating the log-likelihood of the non-gaussian state
  space object. Default is 0, which gives good starting
  value for optimisation. Only used in case of non-Gaussian
  state space model.}

  \item{antithetics}{Logical. If TRUE, two antithetic
  variables are used in simulations, one for location and
  another for scale. Default is TRUE. Only used in case of
  non-Gaussian state space model.}

  \item{taylor}{Logical. If TRUE, control variable based on
  Taylor approximation is used. Default is TRUE. Only used
  in case of non-Gaussian state space model.}

  \item{theta}{Initial values for conditional mode theta.
  Default is \code{object$y}. Only used in case of
  non-Gaussian state space model.}

  \item{maxiter}{Maximum number of iterations used in
  linearisation. Only used in case of non-Gaussian state
  space model.}

  \item{...}{Optional arguments for functions \code{optim}
  and \code{modFun}.}
}
\value{
  A list with elements \item{optim.out}{Output from
  function \code{optim}. } \item{model}{Model with
  estimated parameters. }
}
\description{
  Function \code{fitSSM} finds the maximum likelihood
  estimates for unknown parameters of an arbitary state
  space model if an user defined model building function is
  defined. As a default, \code{fitSSM} estimates the
  non-zero elements, which are marked as NA, of the
  time-invariant covariance matrices H and Q of the given
  model.
}

