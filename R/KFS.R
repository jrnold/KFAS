#' Kalman Filter and Smoother with Exact Diffuse Initialization for Exponential Family State Space Models
#'
#' Performs Kalman filtering and smoothing with exact diffuse initialization
#' using univariate approach for exponential family state space models. 
#' For non-Gaussian models, state smoothing is provided with additional smoothed mean and variance of observations. 
#'
#' Notice that in case of multivariate observations, \code{v}, \code{F}, \code{Finf}, \code{K} and \code{Kinf} are
#' usually not the same as those calculated in usual multivariate Kalman filter. As filtering is done one observation element at the time, 
#' the elements of prediction error \eqn{v_t}{v[t]} are uncorrelated, and \code{F}, \code{Finf}, \code{K} and \code{Kinf} contain only the diagonal elemens of the corresponding covariance matrices.
#'
#' In rare cases of a very long diffuse initialization phase with highly correlated states, cumulative rounding errors
#' in computing \code{Finf} and \code{Pinf} can sometimes cause the diffuse phase end too early.
#' Changing the tolerance parameter \code{tolF} to smaller (or larger) should help.
#'
#' @useDynLib KFAS kfilter statesmooth distsmooth allsmooth
#' @export
#' @param object Object of class \code{SSModel} or \code{KFS} (in which case only smoothing is performed).
#' @param smoothing Perform state or disturbance smoothing or both. Default is \code{"state"} for Gaussian models. For non-Gaussian models, state smoothing is always performed.
#' @param simplify If FALSE, KFS returns some generally not so interesting variables from filtering and smoothing. Default is TRUE.
#' @param transform How to transform the model in case of non-diagonal
#' covariance matrix \eqn{H}. Defaults to \code{"ldl"}. See function \code{\link{transformSSM}} for
#' details.
#' @param nsim Number of independent samples. Default is 100. Only used for non-Gaussian model.
#' @param theta Initial values for conditional mode theta. Default is \code{log(mean(y/u))} for Poisson and 
#' \code{log(mean(y/(u-y)))} for Binomial distribution (or \code{log(mean(y))} in case of \eqn{u_t-y_t = 0}{u[t]-y[t] = 0} for some \eqn{t}). Only used for non-Gaussian model.
#' @param maxiter Maximum number of iterations used in linearisation. Default is 100. Only used for non-Gaussian model.
#' 
#' @return For Gaussian model, a list with the following components: 
#' \item{model}{Original state space model.  }
#' \item{KFS.transform}{Type of H after possible transformation.  }
#' \item{logLik}{Value of the log-likelihood function.  } 
#' \item{a}{One step predictions of states, \eqn{a_t=E(\alpha_t | y_{t-1}, \ldots , y_{1})}{a[t]=E(\alpha[t] | y[t-1], \ldots , y[1])}.  } 
#' \item{P}{Covariance matrices of predicted states, \eqn{P_t=Cov(\alpha_t | y_{t-1}, \ldots , y_{1})}{P[t]=Cov(\alpha[t] | y[t-1], \ldots , y[1])}.  } 
#' \item{Pinf}{Diffuse part of \eqn{P_t}{P[t]}. } 
#' \item{v}{Prediction errors \eqn{v_{i,t} = y_{i,t} - Z_{i,t}a_{i,t}, i=1,\ldots,p}{v[i,t] = y[i,t] - Z[i,t]a[i,t], i=1,\ldots,p},
#' 
#' 
#' where \eqn{a_{i,t}=E(\alpha_t | y_{i-1,t}, \ldots, y_{1,t}, \ldots , y_{1,1})}{a[i,t]=E(\alpha[t] | y[i-1,t], \ldots, y[1,t], \ldots , y[1,1])}.  } 
#' 
#' \item{F}{Prediction error variances \eqn{Var(v_t)}{Var(v[t])}.  } 
#' \item{Finf}{Diffuse part of \eqn{F_t}{F[t]}.  } 
#' \item{d}{The last index of diffuse phase, i.e. the non-diffuse phase began from time \eqn{d+1}.  } 
#' \item{j}{The index of last \eqn{y_{i,t}} of diffuse phase.  } 
#' \item{alphahat}{Smoothed estimates of states, \eqn{E(\alpha_t | y_1, \ldots , y_n)}{E(\alpha[t] | y[1], \ldots , y[n])}. Only computed if \code{smoothing="state"} or \code{smoothing="both"}.  } 
#' \item{V}{Covariances \eqn{Var(\alpha_t | y_1, \ldots , y_n).}{Var(\alpha[t] | y[1], \ldots , y[n]).} Only computed if \code{smoothing="state"} or \code{smoothing="both"}.  } 
#' \item{etahat}{Smoothed disturbance terms \eqn{E(\eta_t | y_1, \ldots , y_n)}{E(\eta[t] | y[1], \ldots , y[n])}.Only computed if \code{smoothing="disturbance"} or \code{smoothing="both"}.  } 
#' \item{V_eta}{Covariances  \eqn{Var(\eta_t | y_1, \ldots , y_n)}{Var(\eta[t] | y[1], \ldots , y[n])}. Only computed if \code{smoothing="disturbance"} or \code{smooth="both"}.  } 
#' \item{epshat}{Smoothed disturbance terms \eqn{E(\epsilon_{t} | y_1, \ldots , y_n)}{E(\epsilon[t] | y[1], \ldots , y[n])}. Only computed if \code{smoothing="disturbance"} or \code{smoothing="both"}. } 
#' \item{V_eps}{Diagonal elements of \eqn{Var(\epsilon_{t} | y_1, \ldots , y_n)}{Var(\epsilon[t] | y[1], \ldots , y[n])}. Note that due to the diagonalization, off-diagonal elements are zero. Only computed if \code{smoothing="disturbance"} or \code{smoothing="both"}.  }
#' In addition, if argument \code{simplify=FALSE}, list contains following components:
#' \item{K}{Covariances \eqn{Cov(\alpha_{t,i}, y_{t,i} | y_{i-1,t}, \ldots, y_{1,t}, y_{t-1}, \ldots , y_{1}), \quad i=1,\ldots,p}{Cov(\alpha[t,i], y[t,i] | y[i-1,t], \ldots, y[1,t], y[t-1], \ldots , y[1]), i=1,\ldots,p}.  }
#' \item{Kinf}{Diffuse part of \eqn{K_t}{K[t]}.  } 
#' \item{r}{Weighted sums of innovations \eqn{v_{t+1}, \ldots , v_{n}}{v[t+1], \ldots , v[n]}.  Notice that in literature t in \eqn{r_t}{r[t]} goes from \eqn{0, \ldots, n}. Here \eqn{t=1, \ldots, n+1}. Same applies to all r and N variables.  } 
#' \item{r0, r1}{Diffuse phase decomposition of \eqn{r_t}{r[t]}.  } 
#' \item{N}{Covariances \eqn{Var(r_t)}{Var(r[t])} .  } 
#' \item{N0, N1, N2}{Diffuse phase decomposition of \eqn{N_t}{N[t]}.   }
#' @return For non-Gaussian model, a list with the following components: 
#' \item{model}{Original state space model with additional elements from function \code{approxSSM}.  }
#' \item{alphahat}{Smoothed estimates of states \eqn{E(\alpha_t | y_1, \ldots , y_n)}{E(\alpha[t] | y[1], \ldots , y[n])}.  } 
#' \item{V}{Covariances  \eqn{Var(\alpha_t | y_1, \ldots , y_n)}{Var(\alpha[t] | y[1], \ldots , y[n])}.  } 
#' \item{yhat}{A time series object containing smoothed means of observation distributions, with parameter \eqn{u_texp(\hat\theta_t)}{u[t]exp(thetahat[t])} for Poisson and \eqn{u_texp(\hat\theta_t)/(1+exp(\hat\theta_t))}{u[t]exp(thetahat[t])/(1+exp(thetahat[t])}.  } 
#' \item{V.yhat}{a vector of length containing smoothed variances of observation distributions.  } 

#' @references Koopman, S.J. and Durbin J. (2000).  Fast filtering and
#' smoothing for non-stationary time series models, Journal of American
#' Statistical Assosiation, 92, 1630-38.  \cr
#'
#' Koopman, S.J. and Durbin J. (2001).  Time Series Analysis by State Space
#' Methods. Oxford: Oxford University Press.  \cr
#'
#' Koopman, S.J. and Durbin J. (2003).  Filtering and smoothing of state vector
#' for diffuse state space models, Journal of Time Series Analysis, Vol. 24,
#' No. 1.  \cr
#'

KFS <- function(object, smoothing = c("state", "disturbance", "both","none"), simplify=TRUE,
        transform = c("ldl", "augment"),nsim=100,theta=NULL,maxiter=100) {
    
    smoothing <- match.arg(arg = smoothing, choices = c("state", "disturbance", "both","none"))
    transform <- match.arg(arg = transform, choices = c("ldl", "augment"))
    
    
    if(!inherits(object,"SSModel") & !inherits(object,"KFS"))
        stop("model must belong to a class 'SSModel' 'KFS'")   
    
    if(inherits(object,"SSModel")){          
        
        if(object$distribution=="Gaussian"){             
                              
            out <- NULL
            out$model <-object
            out$KFS.transform<-"none"
            
            if (object$H_type == "Untransformed") {       
                object <- transformSSM(object, type = transform)
                out$KFS.transform<-object$H_type
            }
            
            tv <- array(0, dim = 5)
            tv[1] <- dim(object$Z)[3] > 1
            tv[2] <- dim(object$H)[3] > 1
            tv[3] <- dim(object$T)[3] > 1
            tv[4] <- dim(object$R)[3] > 1
            tv[5] <- dim(object$Q)[3] > 1    
            
            ymiss <- array(is.na(object$y),dim=c(object$n,object$p))
            storage.mode(ymiss)<-"integer"      
            
            kfout <- .Fortran("kfilter", NAOK = TRUE, array(object$y,dim=c(object$n,object$p)), ymiss, 
                    as.integer(tv), object$Z, object$H, object$T, object$R, object$Q, 
                    object$a1, P1 = object$P1, object$P1inf, object$p, object$n, object$m, 
                    object$k, d = integer(1), j = integer(1), a = array(0, dim = c(object$m, object$n + 1)), 
                    P = array(0, dim = c(object$m, object$m, object$n + 1)), v = array(0, dim = c(object$p, object$n)), 
                    F = array(0, dim = c(object$p, object$n)), K = array(0, dim = c(object$m, object$p, object$n)), 
                    Pinf = array(0, dim = c(object$m, object$m, object$n + 1)), Finf = array(0, dim = c(object$p, object$n)), 
                    Kinf = array(0, dim = c(object$m, object$p, object$n)), lik = double(1), object$tolF, as.integer(sum(object$P1inf)))
            
            
            if (kfout$d == object$n) 
                warning("Degenerate model (d=n), try changing the tolerance parameter tolF or the model.")
            if (kfout$d > 0 & object$m > 1 & !isTRUE(all.equal(min(apply(kfout$Pinf, 3, diag)),0)))
                warning("Possible error in diffuse filtering: Negative variances in Pinf, try changing the tolerance parameter tolF or the model.")
            
            kfout$lik <- kfout$lik - 0.5 * sum(!ymiss) * log(2 * pi) 
            
            ymiss<-t(ymiss)==1
            x<-array(FALSE,c(object$m,object$p,object$n))
            for(i in 1:object$m)
                x[i,,]<-ymiss
            kfout$v[ymiss]<-kfout$F[ymiss]<-kfout$K[x]<-NA
                           
            kfout$Pinf <- kfout$Pinf[1:object$m, 1:object$m, 1:(kfout$d + 1), drop = FALSE]
            if (kfout$d > 0) {
                kfout$Finf[ymiss]<-kfout$Kinf[x]<-NA
                kfout$Finf <- kfout$Finf[, 1:kfout$d, drop = FALSE]
                kfout$Kinf <- kfout$Kinf[, , 1:kfout$d, drop = FALSE]
            } else {
                kfout$Finf <- kfout$Kinf <- NA
            }
            
            rownames(kfout$a)<-rownames(object$a1)
            
            out <- c(out, list(logLik = kfout$lik, a = kfout$a, P = kfout$P, Pinf = kfout$Pinf, 
                            v = kfout$v, F = kfout$F, Finf = kfout$Finf, K = kfout$K, Kinf = kfout$Kinf, 
                            d = kfout$d, j = kfout$j))                                
        } else {
            
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
                   
            tv <- array(0, dim = 5)
            tv[1] <- dim(object$Z)[3] > 1
            tv[2] <- 1
            tv[3] <- dim(object$T)[3] > 1
            tv[4] <- dim(object$R)[3] > 1
            tv[5] <- dim(object$Q)[3] > 1
            
            ymiss <- array(is.na(object$y),dim=c(object$n,object$p))
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
            u <- rnorm(dfu * nsim, mean = 0, sd = 1)
            
            if (dfeps > 0) 
                epsplus[x] <- u[1:(dfeps * nsim)]
            if (dfeta > 0) 
                etaplus[x2] <- u[(dfeps * nsim + 1):(dfeps * nsim + dfeta * nsim)]
            if (nnd > 0) 
                aplus1[nde, ] <- u[(dfeps * nsim + dfeta * nsim + 1):(dfu * nsim)]
            
            
            for (i in 1:nsim) {
                u <- c(etaplus[, , i], epsplus[, , i], aplus1[, i])
                c2[i] <- t(u) %*% c(u)
            }
            q <- pchisq(c2, df = dfu)
            c2 <- sqrt(qchisq(1 - q, dfu)/c2)           
            yo<-which(!is.na(object$y))
            nsim2 <- as.integer(4*nsim)
            
            smoothout <- .Fortran("ngsmooth", PACKAGE = "KFAS", NAOK = TRUE, array(object$y,dim=c(object$n,object$p)), ymiss, 
                    as.integer(tv), object$Z, object$T, object$R, object$Q,  
                    object$a1, object$P1, object$P1inf, object$u, theta=theta,which(object$dist == c("Poisson", "Binomial")),
                    object$p, object$n, object$m, object$k, object$tolF, 
                    as.integer(sum(object$P1inf)), as.integer(nnd), 
                    as.integer(nsim), epsplus, etaplus, aplus1, 
                    c2, object$tol0, info = integer(1), 
                    yo, as.integer(length(yo)), nsim2, 
                    alphahat = array(0,c(object$m,object$n)), V = array(0, c(object$m, object$m, object$n)),
                    meany=array(0,dim=object$n),vary=array(0,dim=object$n),maxiter=as.integer(maxiter),convtol=1e-6,as.integer(nd),as.integer(length(nd)))
            if(maxiter==smoothout$maxiter)
                warning("Maximum number of iterations reached, the linearization did not converge.")
            
            out<-NULL
            out$model <- object
            out$alphahat <- smoothout$alphahat
            rownames(out$alphahat)<-rownames(object$a1)
            out$V <- smoothout$V
            out$yhat<-smoothout$meany
            attributes(out$yhat)<-attributes(object$y)
            out$V.yhat<-smoothout$vary            
            smoothing<-"none" 
        }
    }
    
    if (smoothing != "none") {
        
        if(inherits(object,"KFS")){
            if(object$model$distribution!="Gaussian")
                stop("Input object already contains smoothed states.")
            out<-object
            object<-out$model
            if(out$KFS.transform!="none"){
                object <- transformSSM(object, type = transform)
            }
        }            
        
        tv <- array(0, dim = 5)
        tv[1] <- dim(object$Z)[3] > 1
        tv[2] <- dim(object$H)[3] > 1
        tv[3] <- dim(object$T)[3] > 1
        tv[4] <- dim(object$R)[3] > 1
        tv[5] <- dim(object$Q)[3] > 1    
        
        ymiss <- array(is.na(object$y),dim=c(object$n,object$p))
        storage.mode(ymiss)<-"integer"
        
        if(simplify==FALSE | smoothing=="both"){
            
            smoothout <- .Fortran("allsmooth", NAOK = TRUE, ymiss, 
                    as.integer(tv), object$Z, object$H,object$T, object$R, object$Q,
                    object$p, object$n, object$m, object$k, out$d, out$j,
                    out$a, out$P, out$v, out$F, out$K, r = array(0, dim = c(object$m, object$n + 1)), 
                    r0 = array(0, dim = c(object$m, out$d + 1)), r1 = array(0, dim = c(object$m, out$d + 1)), 
                    N = array(0, dim = c(object$m, object$m, object$n + 1)), 
                    N0 = array(0, dim = c(object$m, object$m, out$d + 1)), 
                    N1 = array(0, dim = c(object$m, object$m, out$d + 1)), 
                    N2 = array(0, dim = c(object$m, object$m, out$d + 1)), out$Pinf, out$Kinf, out$Finf, 
                    object$tolF, alphahat = array(0, dim = c(object$m, object$n)), 
                    V = array(0, dim = c(object$m, object$m, object$n)), 
                    epshat = array(0, dim = c(object$p, object$n)),V_eps = array(0, dim = c(object$p, object$n)),
                    etahat = array(0, dim = c(object$k,object$n)), V_eta = array(0, dim = c(object$k, object$k, object$n)),
                    as.integer(out$KFS.transform == "Augmented"))                       
            
            if(smoothing=="state" | smoothing=="both"){
                out$alphahat <- smoothout$alphahat
                out$V <- smoothout$V
                rownames(out$alphahat)<-rownames(object$a1)
            }
            if(smoothing=="disturbance" | smoothing=="both"){
                out$etahat <- smoothout$etahat
                out$V_eta <- smoothout$V_eta
                if(out$KFS.transform != "Augmented"){
                    out$epshat <- smoothout$epshat
                    out$V_eps <- smoothout$V_eps
                }
            }
            if(simplify==FALSE){
                out$r <- smoothout$r
                out$r0 <- smoothout$r0
                out$r1 <- smoothout$r1
                out$N <- smoothout$N
                out$N0 <- smoothout$N0
                out$N1 <- smoothout$N1
                out$N2 <- smoothout$N2   
            }
            
            
        } else {
            if (smoothing == "state") {            
                smoothout <- .Fortran("statesmooth", NAOK = TRUE, ymiss, 
                        as.integer(tv), object$Z, object$T, object$p, object$n, object$m, out$d, out$j,
                        out$a, out$P, out$v, out$F, out$K, out$Pinf, out$Kinf, out$Finf, object$tolF, 
                        alphahat = array(0, dim = c(object$m, object$n)), V = array(0, dim = c(object$m, object$m, object$n)))
                
                out$alphahat <- smoothout$alphahat
                out$V <- smoothout$V           
                rownames(out$alphahat)<-rownames(object$a1)
            } else {                           
                smoothout <- .Fortran("distsmooth", NAOK = TRUE, ymiss, as.integer(tv), object$Z, object$H, 
                        object$T, object$R, object$Q, object$p, object$n, object$m, object$k,
                        out$d, out$j, out$v, out$F, out$K, out$Kinf, out$Finf, 
                        object$tolF, epshat = array(0, dim = c(object$p, object$n)), 
                        V_eps = array(0, dim = c(object$p, object$n)), 
                        etahat = array(0, dim = c(object$k,object$n)), 
                        V_eta = array(0, dim = c(object$k, object$k, object$n)),as.integer(out$KFS.transform == "augment"))
                
                out$etahat <- smoothout$etahat
                out$V_eta <- smoothout$V_eta
                if(out$KFS.transform != "Augmented"){
                    out$epshat <- smoothout$epshat
                    out$V_eps <- smoothout$V_eps
                }
                
            } 
        }        
        
    }
    
    
    class(out)<-"KFS"
    out
} 
