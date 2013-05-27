#' Linear Gaussian Approximation for Non-Gaussian State Space Model
#'
#' Function \code{approxSMM} computes the linear Gaussian approximation of a
#' state space model where the observations have a non-Gaussian exponential family distribution. 
#' Currently only Poisson and Binomial distributions are supported. 
#' 
#' The linear Gaussian approximating model is a model defined by
#' \deqn{\tilde y_t = Z_t \alpha_t + \epsilon_t, \quad \epsilon_t \sim N(0,\tilde H_t),}{ytilde[t] = Z[t]\alpha[t] + \epsilon[t], \epsilon[t] ~ N(0,Htilde[t]),}
#' \deqn{\alpha_{t+1} = T_t \alpha_t + R_t \eta_t, \quad \eta_t \sim N(0,Q_t),}{\alpha[t+1] = T[t]\alpha[t] + R[t]\eta[t], \eta[t] ~ N(0,Q[t]),} 
#' and \eqn{\alpha_1 \sim N(a_1,P_1)}{\alpha[1] ~ N(a[1],P[1])}, where \eqn{\tilde y}{ytilde} and \eqn{\tilde H}{Htilde} is chosen in a way that the linear 
#' Gaussian approximating model has the same conditional mode of \eqn{\theta=Z\alpha} given the observations \eqn{y} 
#' as the original non-Gaussian model. Models also have same curvature at the mode.
#' 
#' The linearization of the exponential family state space model is based on the first two derivatives of the observational logdensity. 
#'
#' The approximating Gaussian model is used in computation of the log-likelihood of the non-Gaussian model and in importance sampling of non-Gaussian model.  
#' 
#' @seealso Importance sampling of non-Gaussian state space models \code{\link{importanceSSM}}, construct a \code{SSModel} object \code{\link{SSModel}}.
#' @export
#' @param object Non-Gaussian state space model object of class \code{SSModel}.
#' @param theta Initial values for conditional mode theta. Default is \code{log(mean(y/u))} for Poisson and 
#' \code{log(mean(y/(u-y)))} for Binomial distribution (or \code{log(mean(y))} in case of \eqn{u_t-y_t = 0}{u[t]-y[t] = 0} for some \eqn{t}).
#' @param maxiter Maximum number of iterations used in linearisation. Default is 100.
#' @return An object which contains the approximating Gaussian state space model with 
#' additional components \code{original.distribution}, \code{original.y}, \code{thetahat}, and \code{iterations} (the number of iterations used).
approxSSM <- function(object, theta=NULL, maxiter = 100) {
    
  ymiss<-array(is.na(object$y),dim=c(object$n,object$p))
  storage.mode(ymiss)<-"integer"
  
    tv <- array(0, dim = 5)
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
    out <- .Fortran("approx", PACKAGE = "KFAS", NAOK = TRUE, array(object$y,dim=c(object$n,object$p)), ymiss, 
        as.integer(tv), object$Z, object$T, object$R, Htilde = array(0,c(1,1,object$n)), object$Q, object$a1, object$P1, 
        object$P1inf, object$p, object$m, 
        object$k, object$n, theta = theta, object$u, ytilde = array(0, dim = c(object$n,1)), 
        which(object$dist == c("Poisson", "Binomial")),  maxiter=as.integer(maxiter), 
        object$tolF, as.integer(sum(object$P1inf)), convtol = 1e-6)
    
    if(maxiter==out$maxiter)
        warning("Maximum number of iterations reached, the linearization did not converge.")
    object$original.distribution<-object$distribution
    object$original.y<-object$y
    
    object$y <- as.ts(out$ytilde)
    object$y[ymiss] <- NA
    attributes(object$y)<-attributes(object$original.y)
    object$H <- out$Htilde
    object$H_type<-"Diagonal"
    object$thetahat<-out$theta
    object$iterations <- out$maxiter
    
    object$distribution<-"Gaussian"
    class(object)<-c("approxSSM","SSModel")
    invisible(object)
}