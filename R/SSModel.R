#' Create a State Space Model Object of class \code{SSModel}
#'
#' Function \code{SSModel} creates a state space object object of class \code{SSModel}
#' which can be used as an input object for various functions of \code{KFAS} package.
#' 
#' The custom state space model is constructed by using the given system matrices \code{Z}, \code{H}, \code{T}, \code{R}, \code{Q}, 
#' \code{a1}, \code{P1} and \code{P1inf}.
#' Matrix or scalar \code{Z} (array in case of time-varying \code{Z}) is used to determine the number of states \eqn{m}.  
#' If some of the other elements of the object are missing, \code{SSModel} uses default values which are identity matrix for 
#' \code{T}, \code{R} (or \eqn{k} first columns of identity matrix) and \code{P1inf}, and zero matrix for \code{H}, \code{Q}, \code{P1} and , \code{a1}. 
#' If \code{P1} is given and \code{P1inf} is not, the it is assumed to be zero matrix. If \code{Q} is given, 
#' it is used to define \eqn{r}, the dimensions of \code{Q}, which can be smaller than \eqn{m} (defaults to \eqn{m}).
#'
#' The linear Gaussian state space model is given by
#' 
#' \deqn{y_t = Z_t \alpha_t + \epsilon_t,}{y[t] = Z[t]\alpha[t] + \epsilon[t], (observation equation)}
#' 
#' \deqn{\alpha_{t+1} = T_t \alpha_t + R_t \eta_t,}{\alpha[t+1] = T[t]\alpha[t] + R[t]\eta[t], (transition equation)}
#' 
#' where \eqn{\epsilon_t \sim N(0,H_t)}{\epsilon[t] ~ N(0,H[t])}, \eqn{\eta_t \sim N(0,Q_t)}{\eta[t] ~ N(0,Q[t])} 
#' and \eqn{\alpha_1 \sim N(a_1,P_1)}{\alpha[1] ~ N(a[1],P[1])} independently of each other. In case of non-Gaussian observations,
#' the observation equation is of form \eqn{p(y_t|\theta_t) = p(y_t|Z_t\alpha_t)}{p(y[t]|\theta[t]) = p(y[t]|Z[t]\alpha[t])},
#' with \eqn{p(y_t|\theta_t)}{p(y[t]|\theta[t])} being one of the following:
#'
#' If observations are Poisson distributed, parameter of Poisson distribution
#' is \eqn{u_t\lambda_t}{u[t]\lambda[t]} and \eqn{\theta_t = log(\lambda_t)}{\theta[t]=log(\lambda[t])}. 
#'
#' If observations are from binomial distribution, \eqn{u} is a vector
#' specifying number the of trials at times \eqn{1,\ldots,n}, and \eqn{\theta_t =
#' log[\pi_t/(1-\pi_t)]}{\theta[t] = log(\pi[t]/(1-\pi[t]))}, where \eqn{\pi_t}{pi[t]} is the probability of success at time \eqn{t}.
#'
#' For non-Gaussian models \eqn{u_t=1}{u[t]=1} as a default. For Gaussian models, parameter is omitted.
#'
#' Only univariate observations are supported when observation equation is non-Gaussian.
#'
#' @export
#' @seealso \code{\link{arimaSSM}} for state space representation of ARIMA model, \code{\link{regSSM}} 
#' for state space representation of a regression model, and \code{\link{structSSM}} for structural time series model.
#' @param y A time series object of class \code{ts}, or a object that can be coerced to such.
#' @param Z System matrix or array of observation equation.
#' @param H Covariance matrix or array of disturbance terms \eqn{\epsilon_t}{\epsilon[t]} of observation equation. Omitted in case of non-Gaussian distributions. Augment the state vector if you want to add additional noise.
#' @param T System matrix or array of transition equation.
#' @param R System matrix or array of transition equation.
#' @param Q Covariance matrix or array of disturbance terms \eqn{\eta_t}{\eta[t]}.
#' @param a1 Expected value of the initial state vector \eqn{\alpha_1}{\alpha[1]}.
#' @param P1 Covariance matrix of \eqn{\alpha_1}{\alpha[1]}.  In the diffuse case the non-diffuse part of \eqn{P_1}{P[1]}.
#' @param P1inf Diffuse part of \eqn{P_1}{P[1]}. Diagonal matrix with ones on diagonal elements which correspond to the unknown initial states.
#' @param u Only used with non-Gaussian distribution. See details.
#' 
#' @param distribution Specify the distribution of the observations. Default is "Gaussian".
#' @param transform The functions of \code{KFAS} require diagonal covariance matrix \eqn{H_t}{H[t]}. If \eqn{H_t}{H[t]} is not diagonal, model can be transformed using one of the two options.
#' Option \code{"ldl"} performs LDL decomposition for covariance matrix \eqn{H_t}{H[t]}, and multiplies the observation equation with the \eqn{L_t^{-1}}{L[t]^{-1}}, so
#' \eqn{\epsilon_t \sim N(0,D_t)}{\epsilon[t] ~ N(0,D[t])}. Option \code{"augment"} adds \eqn{\epsilon_t}{\epsilon[t]} to the state vector, when
#' \eqn{Q_t}{Q[t]} becomes block diagonal with blocks \eqn{Q_t}{Q[t]} and \eqn{H_t}{H[t]}. 
#' In case of univariate series, option \code{"ldl"} only changes the \code{H_type} argument of the model to \code{"Diagonal"}. Default is \code{"none"} 
#' which does no transformation but checks if \eqn{H} is diagonal. If not, \code{H_type} is set to \code{"Untransformed"}.
#' @param tolF Tolerance parameter for Finf.  Smallest value not counted for zero.
#' @param tol0 Tolerance parameter for LDL decomposition, determines which diagonal values are counted as zero.
#' @return object of class \code{SSModel} with elements 
SSModel <- function(y, Z = NULL, H = NULL, T = NULL, 
        R = NULL, Q = NULL, a1 = NULL, P1 = NULL, P1inf = NULL, u = NULL, 
        distribution = c("Gaussian", "Poisson", "Binomial"), 
        transform = c("none", "ldl", "augment"), tolF = .Machine$double.eps^0.5, tol0 = .Machine$double.eps^0.5) {    
    
    transform <- match.arg(arg=transform, choices = c("none", "ldl", "augment"))    
    distribution <- match.arg(arg=distribution, choices = c("Gaussian", "Poisson", "Binomial"))    
    
    if(!(is.null(Z) | is.array(Z))){
        if(is.vector(Z) & length(Z)>1){
            stop("Error in defining the number of states, Z must be either vector of length 1, data frame, array or NULL.")
        } else {
            Z<-data.matrix(Z)
        }
    }    
    
    y<-as.ts(y)
    if(is.array(y)){
        p <- as.integer(dim(y)[2])
        n <- as.integer(dim(y)[1])  
    } else {
        p <- as.integer(1)
        n <- as.integer(length(y))  
    }     
    storage.mode(y)<-"double"
    #define number of states
    
    m <- as.integer(dim(Z)[2]) #number of custom states
    
    
    Z<-array(Z, dim=c(p, m, (n-1) * (max(dim(Z)[3], 0,na.rm=TRUE) > 1) + 1))
    
    if(is.null(T))  T <- diag(m)
    T<-array(T, dim=c(m, m, (n-1) * (max(dim(T)[3], 0,na.rm=TRUE) > 1) + 1))
    if(distribution!="Gaussian"){
        H_type<-"Omitted"
        H <- array(0, dim=c(p, p, 1))
    } else {
        if(is.null(H)){
            H <- array(0, dim=c(p, p, 1))
            H_type<-"Diagonal"
        } else {            
            H <- array(H, dim=c(p, p,  (n-1) * (max(dim(H)[3], 0,na.rm=TRUE) > 1) + 1))
            if(sum(is.na(H))>0){
                H_type<-"Untransformed"
            } else {
                if(transform=="none"){
                    if(p==1){
                        H_type<-"Diagonal"
                    } else {      
                        H_type<-"Diagonal"
                        for(i in 1:dim(H)[3]){
                            if(max(abs(H[,,i][-which(diag(p)==1)]), na.rm=TRUE)>0){
                                H_type<-"Untransformed"
                                break
                            }
                        }
                    }
                } else{
                    H_type<-NULL
                } 
            }
        }   
    }
    
    if(is.null(Q)){
        r<-as.integer(1)
        Q <- array(0, dim=c(1, 1, 1))        
    } else {       
        r<-as.integer(max(dim(Q)[1],1))
        Q <- array(Q, dim=c(r, r,  (n-1) * (max(dim(Q)[3], 0,na.rm=TRUE) > 1) + 1))        
    }  
    if(is.null(R)){
        R<-diag(m)[,1:r,drop=FALSE]
        dim(R)<-c(m,r,1)
    } else {
        R <- array(R, dim=c(m, r,  (n-1) * (max(dim(R)[3], 0,na.rm=TRUE) > 1) + 1))        
    }
    
    if(is.null(P1) & is.null(P1inf)){
        P1<-matrix(0,m,m)
        P1inf<-diag(1,m)
    } else {
        if(is.null(P1)){
            P1<-matrix(0,m,m)    
        } else P1<-matrix(P1,m,m)
        
        if(is.null(P1inf)){
            P1inf<-matrix(0,m,m)
        } else P1inf<-matrix(P1inf,m,m)
    }
    if(is.null(a1)){
        a1<-matrix(0,m,1)    
    } else a1<-matrix(a1,m,1)
    a1[which(diag(P1inf)>0),1]<-0
    
    if(!identical(distribution,"Gaussian")){
        if(is.null(u))
            u<-rep(1,n)    
        u<-array(u,dim=n)    
    }
    if(is.null(rownames(a1)))
        rownames(a1)<-paste0("state",1:m)
    
    object<-list(y=y,Z=Z,H=H,T=T,R=R,Q=Q,a1=a1,P1=P1,P1inf=P1inf,u=as.double(u),p=p,n=n,m=m,k=r,distribution=distribution,H_type=H_type,tolF=tolF,tol0=tol0)
    class(object) <- c("SSModel","customSSM")
    if (transform %in% c("ldl", "augment") & sum(is.na(H))==0) 
        object <- transformSSM(object, type = transform)
    
    invisible(object)
} 
