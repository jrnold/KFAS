#' Maximum Likelihood Estimation of a State Space Model
#'
#' Function \code{fitSSM} finds the maximum likelihood estimates for 
#' unknown parameters of an arbitary state space model if an user defined model building function is defined. As a default,  
#' \code{fitSSM} estimates the non-zero elements, which are marked as NA, of the time-invariant covariance matrices
#' H and Q of the given model.
#'
#'
#' @export
#' @param inits Initial values for \code{optim}
#' @param model Model object of class \code{SSModel}. if \code{ModFun} is defined, this argument is ignored.
#' @param modFun User defined function which builds the model of class \code{SSModel} given the parameters. If NULL, default estimation procedure is used (See details).
#' @param method The method to be used in \code{optim}. Default is \code{"BFGS"}. 
#' @param nsim Number of independent samples used in estimating the
#' log-likelihood of the non-gaussian state space object. Default is 0, which
#' gives good starting value for optimisation. Only used in case of non-Gaussian state space model.
#' @param antithetics Logical. If TRUE, two antithetic variables are used in
#' simulations, one for location and another for scale. Default is TRUE. Only used in case of non-Gaussian state space model.
#' @param taylor Logical. If TRUE, control variable based on Taylor
#' approximation is used. Default is TRUE. Only used in case of non-Gaussian state space model.
#' @param theta Initial values for conditional mode theta. Default is \code{object$y}. Only used in case of non-Gaussian state space model.
#' @param maxiter Maximum number of iterations used in linearisation. Only used in case of non-Gaussian state space model.
#' @param ... Optional arguments for functions \code{optim} and \code{modFun}.
#' @return A list with elements \item{optim.out}{Output from function \code{optim}. }
#' \item{model}{Model with estimated parameters. }
fitSSM <- function(inits, model=NULL, modFun=NULL, method="BFGS", nsim = 0, antithetics = TRUE, taylor = TRUE, theta=NULL, maxiter = 500, ...) {
    
    out<-NULL
    if(!is.null(modFun)){
        likfn1<-function(inits,...){
            model<-modFun(inits,...)
            -logLik(object = model, nsim = nsim, taylor = taylor,  antithetics = antithetics, maxiter = maxiter, theta=theta)
        }
        out$opt <- optim(par = inits, fn = likfn1, method=method, ...)
        out$model<-modFun(out$opt$par,...)         
        
    } else {                 
        
        if(!is.null(model)){     
            if((!is.null(model$H) && dim(model$H)[3]>1) | (dim(model$Q)[3]>1))
                stop("Time-varying covariance matrices are not allowed during the default estimation. Try to build your own modFun function.")
            if((sum(model$Q!=0,na.rm=TRUE) > 0) & (sum(model$H!=0,na.rm=TRUE) >0))
                stop("Fixed values (excluding zeros) are not allowed in covariance matrices during the default estimation. Try to build your own modFun function.")
            original.model<-model
            if(model$distribution=="Gaussian"){
                if(model$H_type!="Augmented")
                    model<-transformSSM(model,type="augment")
                likfn <- function(pars) {
                    model$Q[lower.tri(model$Q[, , 1])] <- 0
                    if(qdparn> 0)                
                        model$Q[, , 1][qd][qdpar] <- exp(pars[1:qdparn])
                    if (qndparn > 0) 
                        model$Q[upper.tri(model$Q[, , 1])][qndpar] <- pars[(qdparn + 1):(qdparn+qndparn)]
                    model$Q[, , 1] <- crossprod(model$Q[, , 1])   
                    if(any(!is.finite(model$Q))) return(0.5*.Machine$double.xmax)
                    model$P1[(original.model$m+1):(original.model$m+original.model$p),(original.model$m+1):(original.model$m+original.model$p)]<-model$Q[(original.model$k+1):(original.model$k+original.model$p),(original.model$k+1):(original.model$k+original.model$p),1]
                    -logLik(object = model)
                }
                
            }else {
                likfn <- function(pars) {
                    model$Q[lower.tri(model$Q[, , 1])] <- 0
                    if(qdparn> 0)                
                        model$Q[, , 1][qd][qdpar] <- exp(pars[1:qdparn])
                    if (qndparn > 0) 
                        model$Q[upper.tri(model$Q[, , 1])][qndpar] <- pars[(qdparn + 1):(qdparn+qndparn)]
                    model$Q[, , 1] <- crossprod(model$Q[, , 1])     
                    if(any(!is.finite(model$Q))) return(0.5*.Machine$double.xmax)
                    -logLik(object = model, nsim = nsim, antithetics = antithetics, taylor = taylor, theta=theta,
                            maxiter = maxiter)
                }
                
            }
            qd<-1L + 0L:(model$k - 1L) * (model$k + 1L)
            qdpar <- which(is.na(model$Q[, , 1][qd]))
            qdparn <- length(qdpar)
            qndpar <- which(is.na(model$Q[upper.tri(model$Q[, , 1])]))
            qndparn <- length(qndpar)
            model$Q[lower.tri(model$Q[, , 1])] <- 0            
            parn<-qdparn+qndparn
            
            if(parn!=length(inits))
                stop("Number of parameters to estimate is not equal to the number of initial values. ")
            
            
            out$opt <- optim(par = inits, fn = likfn, method=method,...)
            pars<-out$opt$par
            model$Q[lower.tri(model$Q[, , 1])] <- 0
            if(qdparn> 0)                
                model$Q[, , 1][qd][qdpar] <- exp(pars[1:qdparn])
            if (qndparn > 0) 
                model$Q[upper.tri(model$Q[, , 1])][qndpar] <- pars[(qdparn + 1):(qdparn+qndparn)]
            model$Q[, , 1] <- crossprod(model$Q[, , 1])
            
            original.model$Q[,,1]<-model$Q[1:original.model$k,1:original.model$k,1]
            if(model$distribution=="Gaussian")
                original.model$H[,,1]<-model$Q[(original.model$k+1):(original.model$k+original.model$p),(original.model$k+1):(original.model$k+original.model$p),1]        
            
            out$model<-original.model
            
        } else stop("Either model or model building function modFun must be provided.")
    }
    out
}
