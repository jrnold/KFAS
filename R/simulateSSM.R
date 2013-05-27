#' Simulation of a Gaussian State Space Model
#'
#' Function \code{simulateSMM} simulates states, disturbances or missing observations of
#' the Gaussian state space object conditionally on the data.
#'
#'
#' Simulation smoother algorithm is from article by J. Durbin and S.J. Koopman
#' (2002).
#'
#' Function can use two antithetic variables, one for location and other for
#' scale, so output contains four blocks of simulated values which correlate
#' which each other (ith block correlates negatively with (i+1)th block, and
#' positively with (i+2)th block etc.).
#'
#' @export
#' @param object Gaussian state space object.
#' @param sim What to simulate. Note that all the simulations are done independently.
#' @param nsim Number of independent samples. Default is 1.
#' @param antithetics Use antithetic variables in simulation. Default is FALSE.
#' @param conditional Simulations are conditional to data. 
#' If FALSE, the initial state \eqn{\alpha_1}{\alpha[1]} is set to \eqn{\hat \alpha_1}{alphahat[1]} computed by \code{KFS}, 
#' and all the observations are removed from the model. Default is TRUE.
#' @references Durbin J. and Koopman, S.J. (2002). A simple and efficient
#' simulation smoother for state space time series analysis, Biometrika, Volume
#' 89, Issue 3
simulateSSM <- function(object, sim = c("states", "disturbances", "observations"), 
    nsim = 1, antithetics = FALSE, conditional=TRUE) {
    sim.what <- match.arg(arg = sim, choices = c("states", "disturbances", "observations"), 
        several.ok = TRUE)
    sims <- NULL
    if(object$distribution!="Gaussian")
        stop("Function is only for Gaussian models.")
    if(!conditional){
        out<-KFS(object,smoothing="state")
        object$y[]<-NA
        object$a1[]<-out$alphahat[,1]
        object$P1inf[]<-0
        object$P1[]<-0
        }
    if ("states" %in% sim.what) 
        sims$states <- simState(object, nsim = nsim, antithetics = antithetics)
    
    if ("disturbances" %in% sim.what) 
        sims$disturbances <- simDisturbance(object, nsim = nsim, antithetics = antithetics)
    
    if ("observations" %in% sim.what) 
        sims$observations <- simObs(object, nsim = nsim, antithetics = antithetics)
    sims
} 
