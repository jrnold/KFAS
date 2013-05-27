#' Estimate the Missing Observations of a State Space Model
#'
#' Function \code{predict.KFS} computes the expected values of missing
#' observations given the observed data.
#' 
#' 
#'
#' @export
#' @useDynLib KFAS predobs
#' @method predict KFS
#' @S3method predict KFS
#' @aliases predict predict.KFS
#' @param object object of class \code{KFS} ie. the output from function KFS.
#' @param fill If FALSE, only predictions of missing observations are returned, and other time points are markes as NA. 
#' This is convinient for plotting purposes. If TRUE, original time series is filled with predicted values for missing observations, 
#' and \eqn{F_t = 0}{F[t] = 0} if observation is not missing. Default is FALSE.
#' @param \dots Ignored.
#' @return A list with the following components: 
#' \item{y}{Time series object with missing observations replaced by \eqn{E(y_t|y)}{E(y[t]|y)}.  } 
#' \item{F}{Covariances \eqn{Cov(y_t|y)}{Cov(y[t]|y)}. 
#' Note that this is the usual multivariate version of \eqn{F_t}{F[t]} given by \eqn{Z_tP_tZ'_t + H_t}{Z[t]P[t]Z'[t] + H[t]}, not the univariate version given by KFS.  }
predict.KFS <- function(object, fill=FALSE, ...) {
    ymiss<-is.na(object$model$y)
    if(sum(ymiss)==0) stop("No missing observations, nothing to predict.")
    
    if(object$model$distribution=="Gaussian"){   
        
        if(is.null(object$alphahat)){
            object<-KFS(object,smoothing="state")
        }
        a<-object$alphahat
        P<-object$V
        
        
        storage.mode(ymiss)<-"integer"
        pred <- .Fortran("predobs", NAOK = TRUE, ymiss, as.integer(dim(object$model$Z)[3] > 1),
                as.integer(dim(object$model$H)[3] > 1),
                object$model$Z,object$model$H,a,P,object$model$p, object$model$n,
                object$model$m,theta=array(object$model$y,c(object$model$n,object$model$p)),
                sigma=array(-1000,c(object$model$p,object$model$p,object$model$n)))
        if(fill){
            pred$sigma[pred$sigma < -100]<-0 
        } else {
            pred$sigma[pred$sigma < -100]<-NA   
            pred$theta[!ymiss]<-NA
        }
        
    } else{
        pred<-NULL
        pred$theta<-object$yhat
        pred$sigma<-object$V.yhat
        if(fill){
            pred$sigma[!ymiss]<-0 
            pred$theta[!ymiss]<-object$model$y[!ymiss]
        } else {
            pred$sigma[!ymiss]<-NA
            pred$theta[!ymiss]<-NA
        }
    }
    attributes(pred$theta)<-attributes(object$model$y)    
    list(y = pred$theta, F = pred$sigma)
}


#' Extract the filtered or smoothed signal of a State Space Model
#'
#' Function \code{signal} extracts the filtered or smoothed signal of a State Space model depending on the input object. 
#'
#' @export
#' @useDynLib KFAS signaltheta signalthetang
#' @param object Object of class \code{KFS}.
#' @param states Which states are combined? Default is NULL which combines all states according to \eqn{Z_t}{Z[t]}.
#' @return
#'\item{signal}{Time series object of filtered signal \eqn{Z_ta_t}{Z[t]a[t]} or smoothed signal \eqn{Z_t\hat\alpha_t}{Z[t]\alpha[t]} using only the defined states. Filtered signal is computed only for non-diffuse phase.  }
#' \item{variance}{Cov(\eqn{Z_ta_t}{Z[t]a[t]}) or Cov(\eqn{Z_t\hat\alpha_t}{Z[t]\alpha[t]}) using only the defined states. Covariances of filtered signal is computed only for non-diffuse phase.  }
signal <- function(object,states=NULL) {
    
    if(!is.null(states)){
        if(!is.numeric(states)){
            which(states==object$model$a1)
        } else states<-as.integer(states)
        if(min(states)<1 | max(states)>object$model$m)
            stop("Vector states should contain the indices of the states which are combined.")
    } else states<-as.integer(1:object$model$m)
    
    #if(object$model$distribution=="Gaussian"){     
    
    if(is.null(object$alphahat)){
        a<-object$a
        P<-object$P
        d<-as.integer(object$d)
    } else {
        a<-object$alphahat
        P<-object$V
        d<-as.integer(0)
    }                   
    
    signal <- .Fortran("signaltheta", NAOK = TRUE, as.integer(dim(object$model$Z)[3] > 1),
            object$model$Z,a[1:object$model$m,1:object$model$n],P[1:object$model$m,1:object$model$m,1:object$model$n],object$model$p, object$model$n,
            object$model$m,theta=array(0,c(object$model$n,object$model$p)),
            sigma=array(0,c(object$model$p,object$model$p,object$model$n)),d=d,states,as.integer(length(states)))
    if(d>0){
        signal$theta[1:d,]<-NA
        signal$sigma[,,1:d]<-NA
    }
    
#    } else{
#        iSample<-importanceSSM(object$model, nsim = nsim, ...)        
#       
#            signal <- .Fortran("signalthetang", NAOK = TRUE, as.integer(dim(object$model$Z)[3] > 1),
#                    object$model$Z,iSample$states,iSample$w/sum(iSample$w), object$model$n,
#                    object$model$m,theta=array(0,c(object$model$n)),
#                    sigma=array(0,c(object$model$n)),states,
#                    as.integer(length(states)),as.integer(length(iSample$w)))
#        }
    
    attributes(signal$theta)<-attributes(object$model$y)   
    list(signal = signal$theta, variance = signal$sigma)
}