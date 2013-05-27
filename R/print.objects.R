#' Print SSModel Object
#' @S3method print SSModel
#' @method print SSModel
#' @param x SSModel object
#' @param ... Ignored.

print.SSModel<-function(x, ...){
    cat("State space model object of class SSModel\n")    
    cat("y = \n")
    if(x$n<10){       
        print(x$y)
    }else {
        cat("First 10 time points:\n")        
        print(window(x$y,start=start(x$y)[1],end=start(x$y)[1]+9))
       
    }
    cat("\nZ = \n")
    if(dim(x$Z)[3]==1){       
        print(matrix(x$Z[,,1],x$p,x$m))
    }else {
        if(x$n<10){
            print(x$Z)
        } else {
            cat("First 10 time points:\n")
            print(x$Z[,,1:10,drop=FALSE])
        }
    }
    cat("\nH = \n")
    if(dim(x$H)[3]==1){       
        print(x$H[,,1])
    }else {
        if(x$n<10){
            print(x$H)
        } else {
            cat("First 10 time points:\n")
            print(x$H[,,1:10,drop=FALSE])
        }
    }    
    
    cat("\nT = \n")
    if(dim(x$T)[3]==1){       
        print(x$T[,,1])
    }else {
        if(x$n<10){
            print(x$T)
        } else {
            cat("First 10 time points:\n")
            print(x$T[,,1:10,drop=FALSE])
        }
    }  
    cat("\nR = \n")
    if(dim(x$R)[3]==1){       
        print(matrix(x$R[,,1],x$m,x$k))
    }else {
        if(x$n<10){
            print(x$R)
        } else {
            cat("First 10 time points:\n")
            print(x$R[,,1:10,drop=FALSE])
        }
    }  
    cat("\nQ = \n")
    if(dim(x$Q)[3]==1){       
        print(x$Q[,,1])
    }else {
        if(x$n<10){
            print(x$Q)
        } else {
            cat("First 10 time points:\n")
            print(x$Q[,,1:10,drop=FALSE])
        }
    }  
    
    cat("\na1 = \n")   
    print(x$a1[,,drop=(x$m==1)])
    
    cat("\nP1 = \n")          
    print(x$P1[,,drop=TRUE])
    
    cat("\nP1inf = \n")       
    print(x$P1inf[,,drop=TRUE])
    
    if(!identical(x$distribution,"Gaussian")){
        cat("\nu = \n")          
        print(x$u)    
    }
}   
#' Print Ouput of Kalman Filter and Smoother
#' @S3method print KFS
#' @method print KFS
#' @param x output object from function KFS.
#' @param ... Ignored.

print.KFS<-function(x, ...){
    cat("_Output of Kalman filtering_\n")    
    cat("\nLog-likelihood =",x$logLik,"\n")          
    
    cat("\n Filtered values of the initial and last states:\n")
    cat("\na_1 =\n")
    print(x$a[,1])  
     
    cat("\na_",x$model$n," =\n",sep="")         
    print(x$a[,x$model$n+1]) 
    
    cat("\nCovariances of the filtered values of initial and last states:\n")
    cat("\nP_1 =\n")  
    print(x$P[,,1])
    cat("\nP_",x$model$n," =\n",sep="")         
    print(x$P[,,x$model$n+1]) 
    
    if(!is.null(x$alpha)){
        cat("\n_Output of state smoothing_\n")
        
        cat("\nSmoothed values of the initial and last states:\n")
        cat("\nalpha_1 =\n")
        print(x$alpha[,1])  
        
        cat("\nalpha_",x$model$n," =\n",sep="")         
        print(x$alpha[,x$model$n]) 
        
        cat("\nCovariances of the smoothed values of initial and last states:\n")
        cat("\nV_1 =\n")  
        print(x$V[,,1])
        cat("\nV_",x$model$n," =\n",sep="")         
        print(x$V[,,x$model$n])        
    }
    if(!is.null(x$epshat)){
        cat("\n_Output of disturbance smoothing_\n")
        
        cat("\nSmoothed values of the initial and last epsilon disturbances:\n")
        cat("\nepshat_1 =\n")
        print(x$epshat[,1])  
        
        cat("\nepshat_",x$model$n," =\n",sep="")         
        print(x$epshat[,x$model$n]) 
        
        cat("\nCovariances of the smoothed values of initial and last epsilon disturbances:\n")
        cat("\nV_eps_1 =\n")  
        print(x$V_eps[,1])
        cat("\nV_epsr_",x$model$n," =\n",sep="")         
        print(x$V_eps[,x$model$n])        
    }
    if(!is.null(x$etahat)){
        cat("\n_Output of disturbance smoothing_\n")
        
        cat("\nSmoothed values of the initial and last eta disturbances:\n")
        cat("\netahat_1 =\n")
        print(x$etahat[,1])  
        
        cat("\netahat_",x$model$n," =\n",sep="")         
        print(x$etahat[,x$model$n]) 
        
        cat("\nCovariances of the smoothed values of initial and last eta disturbances:\n")
        cat("\nV_eta_1 =\n")  
        print(x$V_eta[,,1])
        cat("\nV_eta_",x$model$n," =\n",sep="")         
        print(x$V_eta[,,x$model$n])        
    }

}
#' Extract Standardized Residuals of Kalman Filter and Smoother output
#' @S3method residuals KFS
#' @method residuals KFS
#' @details For object of class KFS, three types of residuals can be computed:
#' Standardized one-step ahead prediction residuals are defined as
#' \deqn{v_{i,t}/\sqrt{F_{i,t}}, \quad i=1,\ldots,p,t=d+1,\ldots,n,}{v[i,t]F[i,t]^-0.5, i=1,\ldots,p, t=d+1,\ldots,n,} 
#' with residuals being undefined in diffuse phase.
#' 
#' Residuals based on the smoothed disturbance terms are defined as
#' \deqn{\hat \epsilon_{i,t}/\sqrt{Var(\epsilon_{i,t})}, \quad i=1,\ldots,p, t=1,\ldots,n,}{\hat \epsilon[i,t]Var(\epsilon[i,t])^-0.5, i=1,\ldots,p, t=1,\ldots,n,} and
#' \deqn{L^{-1}_t \hat \eta_t, \quad t=1,\ldots,n,}{L^{-1}[t] \eta[t], t=1,\ldots,n,} where \eqn{L_t}{L[t]} is the lower triangular matrix from Cholesky decomposition of \eqn{V_{\eta,t}}{V[\eta,t]}
#' @param object KFS object
#' @param ... Ignored.

residuals.KFS<-function(object,...){
    out<-NULL
    out$v<-t(object$v/sqrt(object$F))
    out$v[1:object$d,]<-NA
    attributes(out$v)<-attributes(object$model$y)
    
    if(!is.null(object$epshat)){
        out$eps<-t(object$epshat/sqrt(object$V_eps))        
        attributes(out$eps)<-attributes(object$model$y)
    }
    if(!is.null(object$etahat)){
        eta<-array(0,c(object$model$n,length(object$etahat[,1])))
        for(i in 1:object$model$n){
            x<-try(chol(solve(object$V_eta[,,i]))%*%object$etahat[,i],TRUE)
            if(inherits(x,"try-error")){
                warning(paste("Could not compute the smoothed eta residuals, V_eta[,,",i,"] is not invertible",sep=""))
                break
            } else  eta[i,]<-x          
        } 
        out$eta<-ts(eta)
        tsp(out$eta)<-tsp(out$v)
    }      
    out
    
}