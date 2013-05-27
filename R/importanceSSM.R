#' Importance Sampling of Non-Gaussian State Space Model
#'
#' Importance Sampling of Non-Gaussian State Space Model.
#'
#' Function \code{importanceSSM} simulates states of the non-Gaussian state space model conditioned with the observations,
#' returning the simulated samples of the  states with the importance weights.
#'
#' @export
#' @param model Non-Gaussian state space model object of class \code{SSModel}.
#' @param nsim Number of independent samples. Default is 1000.
#' @param save.model Return the original model with the samples. Default is FALSE.
#' @param theta Initial values for conditional mode theta. Default is \code{log(mean(y/u))} for Poisson and 
#' \code{log(mean(y/(u-y)))} for Binomial distribution (or \code{log(mean(y))} in case of \eqn{u_t-y_t = 0}{u[t]-y[t] = 0} for some \eqn{t}).
#' @param antithetics Logical. If TRUE, two antithetic variables are used in
#' simulations, one for location and another for scale. Default is TRUE.
#' @param maxiter Maximum number of iterations used in linearisation. Default is 100.
importanceSSM <- function(model, nsim = 1000, save.model=FALSE, theta=NULL,
        antithetics = TRUE, maxiter = 100) {

    
    if(is.null(theta)){
        if(model$distribution=="Poisson"){
            theta<-array(log(mean(model$y/model$u,na.rm=TRUE)),dim=model$n)    
        } else {
            if(min(model$u-model$y,na.rm=TRUE)==0){
                theta<-array(log(mean(model$y,na.rm=TRUE)),dim=model$n)
            } else {
                theta<-array(log(mean(model$y/(model$u-model$y),na.rm=TRUE)),dim=model$n)
            }
        }
    } else {
        theta<-array(theta,dim=model$n)
    }
    
    tv <- array(0, dim = 5)
    tv[1] <- dim(model$Z)[3] > 1
    tv[2] <- 1
    tv[3] <- dim(model$T)[3] > 1
    tv[4] <- dim(model$R)[3] > 1
    tv[5] <- dim(model$Q)[3] > 1
    
    ymiss <- array(is.na(model$y),dim=c(model$n,model$p))
    storage.mode(ymiss)<-"integer"
    
    epsplus <- array(0, c(1, model$n, nsim))
    etaplus <- array(0, c(model$k, model$n, nsim))
    aplus1 <- array(0, dim = c(model$m, nsim))
    c2 <- numeric(nsim)
    
    x <- c(!ymiss)
    x <- array(x, c(1, model$n, nsim))
    dfeps <- sum(x)/nsim
    
    x2 <- array(apply(model$Q, 3, diag) > model$tol0, c(model$k, (model$n - 1) * tv[5] + 1))
    x2 <- array(x2, c(model$k, model$n, nsim))
    dfeta <- sum(x2)/nsim
    
    nde <- which(diag(model$P1) > model$tol0)
    nnd <- length(nde)
    nd <- which(diag(model$P1inf) == 0)
    dfu <- dfeps + dfeta + nnd   
    u <- rnorm(dfu * nsim, mean = 0, sd = 1)
    
    if (dfeps > 0) 
        epsplus[x] <- u[1:(dfeps * nsim)]
    if (dfeta > 0) 
        etaplus[x2] <- u[(dfeps * nsim + 1):(dfeps * nsim + dfeta * nsim)]
    if (nnd > 0) 
        aplus1[nde, ] <- u[(dfeps * nsim + dfeta * nsim + 1):(dfu * nsim)]
    
    if (antithetics) {
        for (i in 1:nsim) {
            u <- c(etaplus[, , i], epsplus[, , i], aplus1[, i])
            c2[i] <- t(u) %*% c(u)
        }
        q <- pchisq(c2, df = dfu)
        c2 <- sqrt(qchisq(1 - q, dfu)/c2)
    }    
    
    nsim2 <- as.integer(3 * antithetics * nsim + nsim)
    yo <- which(!is.na(model$y))
    out <- .Fortran("importance", PACKAGE = "KFAS", NAOK = TRUE, array(model$y,dim=c(model$n,model$p)), ymiss, 
            as.integer(tv), model$Z,  model$T, model$R, model$Q, model$a1, model$P1, model$P1inf, model$u, 
            dist = which(model$dist == c("Poisson", "Binomial")), model$p, model$n, model$m, model$k, theta, maxiter=as.integer(maxiter), model$tolF, 
            as.integer(sum(model$P1inf)), convtol = 1e-6, as.integer(nnd), 
            as.integer(nsim), epsplus, etaplus, aplus1, c2, model$tol0, info = integer(1), as.integer(antithetics), 
            yo, as.integer(length(yo)), nsim2, w = numeric(nsim2), asim = array(0, c(model$m, model$n, nsim2)),as.integer(nd),as.integer(length(nd)))
    if(maxiter==out$maxiter)
        warning("Maximum number of iterations reached, the linearization did not converge.")
    
    out <- list(states = out$asim, weights = out$w)     
    if(save.model)
        out$model<-model
    
    class(out)<-"importanceSSM"
    out
}
