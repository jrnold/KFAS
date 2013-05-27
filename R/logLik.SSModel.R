#' Log-likelihood of the State Space Model.
#'
#' Function \code{logLik.SSmodel} computes the log-likelihood value of a state-space
#' model.
#'
#'
#' @useDynLib KFAS glogliku gloglikd gloglika ngloglik
#' @inheritParams approxSSM
#' @export
#' @S3method logLik SSModel
#' @method logLik SSModel
#' @aliases logLik logLik.SSModel
#' @param object State space model of class \code{SSModel}.
#' @param nsim Number of independent samples used in estimating the
#' log-likelihood of the non-gaussian state space model. Default is 0, which
#' gives good starting value for optimization. Only used for non-Gaussian model.
#' @param antithetics Logical. If TRUE, two antithetic variables are used in
#' simulations, one for location and another for scale. Default is TRUE. Only used for non-Gaussian model.
#' @param taylor Logical. If TRUE, control variable based on Taylor
#' series is used. Default is TRUE. Only used for non-Gaussian model.
#' @param theta Initial values for conditional mode theta. Default is \code{log(mean(y/u))} for Poisson and 
#' \code{log(mean(y/(u-y)))} for Binomial distribution (or \code{log(mean(y))} in case of \eqn{u_t-y_t = 0}{u[t]-y[t] = 0} for some \eqn{t}). Only used for non-Gaussian model.
#' @param fix.seed Use fixed seed. If FALSE, no fixed seed is used. If fix.seed
#' is positive value, the value is used as a seed via set.seed function. 
#' Default is TRUE, so that the variation in random number generation does not affect numerical optimization algorithms. Only used for non-Gaussian model.
#' @param ... Ignored.
#' @return \item{}{log-likelihood of the state space model.}
logLik.SSModel <- function(object, nsim = 0, antithetics = TRUE, taylor = TRUE, theta=NULL, 
        maxiter = 100, fix.seed = TRUE, ...) {
    
    ymiss <- array(is.na(object$y),dim=c(object$n,object$p))
    storage.mode(ymiss)<-"integer" 
    if (object$distribution == "Gaussian") {
        
        lik <- 0
        kfout <- NULL
        if (object$p == 1) {
            tv <- array(0, dim = 5)
            tv[1] <- dim(object$Z)[3] > 1
            tv[2] <- dim(object$H)[3] > 1
            tv[3] <- dim(object$T)[3] > 1
            tv[4] <- dim(object$R)[3] > 1
            tv[5] <- dim(object$Q)[3] > 1           
            kfout <- .Fortran("glogliku", PACKAGE = "KFAS", NAOK = TRUE, array(object$y,dim=c(object$n,object$p)), 
                    ymiss, as.integer(tv), object$Z, object$H, object$T, object$R, 
                    object$Q, object$a1, object$P1, object$P1inf,object$m, object$k, 
                    object$n, lik = lik, object$tolF, as.integer(sum(object$P1inf)))
            
        } else {
            if (object$H_type == "Untransformed") 
                object <- transformSSM(object, type = "ldl")
            
            
            tv <- array(0, dim = 5)
            tv[1] <- dim(object$Z)[3] > 1
            tv[2] <- dim(object$H)[3] > 1
            tv[3] <- dim(object$T)[3] > 1
            tv[4] <- dim(object$R)[3] > 1
            tv[5] <- dim(object$Q)[3] > 1
            if (object$H_type != "Augmented") {
                kfout <- .Fortran("gloglikd", PACKAGE = "KFAS", NAOK = TRUE, array(object$y,dim=c(object$n,object$p)), 
                        ymiss, as.integer(tv), object$Z, object$H, object$T, object$R, 
                        object$Q, object$a1, object$P1, object$P1inf, object$p, 
                        object$m, object$k, object$n, lik = lik, object$tolF, as.integer(sum(object$P1inf)))
            } else {
                #augment
                kfout <- .Fortran("gloglika", PACKAGE = "KFAS", NAOK = TRUE, array(object$y,dim=c(object$n,object$p)), 
                        ymiss, as.integer(tv), object$Z, object$T, object$R, object$Q, 
                        object$a1, object$P1, object$P1inf, object$p, object$m, object$k, 
                        object$n, lik = lik, object$tolF, as.integer(sum(object$P1inf)))
            }
        }
        logLik <- kfout$lik - 0.5 * sum(!ymiss) * log(2 * pi)
    } else {
        
        
        
        tv <- array(0, dim = 4)
        tv[1] <- dim(object$Z)[3] > 1
        tv[2] <- 1
        tv[3] <- dim(object$T)[3] > 1
        tv[4] <- dim(object$R)[3] > 1
        tv[5] <- dim(object$Q)[3] > 1
        
        if(is.null(theta)){
            if(object$distribution=="Poisson"){
                theta<-array(log(mean(object$y/object$u,na.rm=TRUE)),dim=object$n)    
            } else {
                if(min(object$u-object$y,na.rm=TRUE)==0){
                    theta<-array(log(mean(object$y,na.rm=TRUE)),dim=object$n)
                } else {
                    theta<-array(log(mean(object$y/(object$u-object$y),na.rm=TRUE)),dim=object$n)
                }
            }
        } else {
            theta<-array(theta,dim=object$n)
        }
        
        if (nsim == 0) {
            nsim <- 1
            sim <- 0
            epsplus <- array(0, c(1, object$n, nsim))
            etaplus <- array(0, c(object$k, object$n, nsim))
            aplus1 <- array(0, dim = c(object$m, nsim))
            c2 <- numeric(nsim)
            nnd <- 0
            nd <- which(diag(object$P1inf) == 0)
        } else {
            sim <- 1
            
            
            ymiss <- is.na(object$y)
            storage.mode(ymiss)<-"integer"
            
            
            epsplus <- array(0, c(1, object$n, nsim))
            etaplus <- array(0, c(object$k, object$n, nsim))
            aplus1 <- array(0, dim = c(object$m, nsim))
            c2 <- numeric(nsim)
            
            x <- c(!ymiss)
            x <- array(x, c(1, object$n, nsim))
            dfeps <- sum(x)/nsim
            
            x2 <- array(apply(object$Q, 3, diag) > object$tol0, c(object$k, (object$n - 1) * tv[5] + 1))
            x2 <- array(x2, c(object$k, object$n, nsim))
            dfeta <- sum(x2)/nsim
            
            nde <- which(diag(object$P1) > object$tol0)
            nnd <- length(nde)
            nd <- which(diag(object$P1inf) == 0)
            dfu <- dfeps + dfeta + nnd
            if (fix.seed > 0) 
                set.seed(fix.seed)
            u <- rnorm(n= dfu * nsim, mean = 0, sd = 1)
            
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
        }
        ##
        nsim2 <- as.integer(max(sim * (3 * antithetics * nsim + nsim), 1))
        yo <- which(!is.na(object$y))
        out <- .Fortran("ngloglik", PACKAGE = "KFAS", NAOK = TRUE, yt=array(object$y,dim=c(object$n,object$p)), 
                ymiss, as.integer(tv), object$Z, Htilde=array(0,c(1,1,object$n)),
                object$T, object$R, object$Q, object$a1, object$P1, object$P1inf, 
                object$p, object$m, object$k, object$n, lik = double(1), theta = theta, 
                object$u, ytilde = array(0, dim = c(object$n, 1)), dist=which(object$distribution == c("Poisson", "Binomial")), as.integer(maxiter), 
                object$tolF, as.integer(sum(object$P1inf)), convtol = 1e-6, taylor=as.integer(taylor), 
                as.integer(nnd), as.integer(nsim), epsplus, etaplus, aplus1, 
                c2, object$tol0, info = integer(1), as.integer(antithetics), as.integer(sim), 
                yo, as.integer(length(yo)), nsim2, eg = double(1),as.integer(nd),as.integer(length(nd)))
        
        if(out$taylor== -1)
            warning("Control variable did not work properly and therefore was not used in computing the log-likelihood.")        
        if (out$dist == 1) 
            constants <- -log(nsim2) + sum(dpois(x = out$yt[yo,1], lambda = object$u * exp(out$theta[yo]), 
                                    log = TRUE)) - sum(dnorm(x = out$ytilde[yo,1], mean = out$theta[yo], 
                                    sd = sqrt(out$Htilde[1, 1, yo]), log = TRUE))
        if (out$dist == 2) 
            constants <- -log(nsim2) + sum(dbinom(x = out$yt[yo,1], size = object$u[yo], 
                                    prob = (exp(out$theta)/(1 + exp(out$theta)))[yo], log = TRUE)) - 
                    sum(dnorm(x = out$ytilde[yo,1], mean = out$theta[yo], sd = sqrt(out$Htilde[1, 1, yo]), log = TRUE))
        
        logLik <- out$lik - 0.5 * sum(!is.na(object$y)) * log(2 * pi) + out$eg + constants  #+sum(log(sqrt(2*pi*out$Ht[1,1,yo])))
        
    }
    logLik
} 
