#' Create a State Space Model Representation of Linear Regression Model
#'
#' Function regSSM creates a state space representation of linear regression model.
#' 
#' The linear Gaussian state space model is given by
#' 
#' \deqn{y_t = X_t \beta_t + \epsilon_t,}{y[t] = Z[t]\alpha[t] + \epsilon[t], (observation equation)}
#' 
#' \deqn{\alpha_{t+1} = T_t \alpha_t + R_t \eta_t,}{\alpha[t+1] = T[t]\alpha[t] + R[t]\eta[t], (transition equation)}
#' 
#' where \eqn{\epsilon_t ~ N(0,H_t)}{\epsilon[t] ~ N(0,H[t])}, \eqn{\eta_t ~ N(0,Q_t)}{\eta[t] ~ N(0,Q[t])} 
#' and \eqn{\alpha_1 ~ N(a_1,P_1)}{\alpha[1] ~ N(a[1],P[1])} independently of each other. In case of non-Gaussian observations,
#' the observation equation is of form \eqn{p(y_t|\theta_t) = p(y_t|Z_t\alpha_t)}{p(y[t]|\theta[t]) = p(y[t]|Z[t]\alpha[t])},
#' with \eqn{p(y_t|\theta_t)}{p(y[t]|\theta[t])} being one of the following:
#'

#'
#' @export
#'  @inheritParams SSModel
#'  @seealso \code{\link{arimaSSM}} for state space representation of ARIMA model, \code{\link{structSSM}} for structural time series model, and \code{\link{SSModel}} for custom \code{SSModel} object.
#' @param X A \eqn{n \times k}{n*k} matrix of explanatory variables, with each column containing one explanatory variable, or a list of length \eqn{p} 
#' containing \eqn{X} matrices for each series. If X is matrix, it is assumed that all \eqn{p} series use same explanatory variables.  
#' @param H A \eqn{p \times p}{p*p} covariance matrix (or \eqn{p \times p \times n}{p*p*n} array in of time-varying case) of the disturbance terms 
#'  \eqn{\epsilon_t}{epsilon[t]} of the observation equation. Default gives \eqn{p \times p}{p*p} zero matrix.
#'  Omitted in case of non-Gaussian distributions. Augment the state vector if you want to add additional noise.
#' @param Q A \eqn{r \times r}{r*r} (or \eqn{r \times r \times n}{r*r*n} array in of time-varying case) covariance matrix of the disturbance terms 
#'  \eqn{\eta_t}{\eta[t]} of the system equation.  Default is \eqn{m \times m}{m*m} zero matrix ie. ordinary time-invariant regression.
regSSM <- function(y, X, H=NULL, Q=NULL, u=NULL, distribution = c("Gaussian", "Poisson", "Binomial"), 
        transform = c("none", "ldl", "augment"), tolF = .Machine$double.eps^0.5, tol0 = .Machine$double.eps^0.5) {
    
    transform <- match.arg(arg=transform, choices = c("none", "ldl", "augment"))    
    distribution <- match.arg(arg=distribution, choices = c("Gaussian", "Poisson", "Binomial"))    
    
    y<-as.ts(y)
    if(is.array(y)){
        p <- as.integer(dim(y)[2])
        n <- as.integer(dim(y)[1])  
    } else {
        p <- as.integer(1)
        n <- as.integer(length(y))  
    }     
    storage.mode(y)<-"double"
    
    if(p>1 & !identical(distribution,"Gaussian")) stop("Only univariate series are supported for non-Gaussian models.")
    #define number of states    
    
    if(!is.list(X)){       
        
        if(!is.ts(X) & (!identical(dim(X)[1],n) | is.null(dim(X)[2]))) stop("X is not an n times k matrix or list of such matrices.")
        X2<-vector("list",p)
        for(i in 1:p)
            X2[[i]]<-data.matrix(X)
        X<-X2
    } 
    xdim <- sapply(X,dim) #number of regressors including the constant
    if(!identical(sum(xdim[1,]==n),p)) stop("X is not an n times k matrix or list of such matrices.")
    
    kn <- xdim[2,] #number of regressors including the constant
    m <- as.integer(sum(kn)) # number of states
    tspy<-attributes(y)
    y<-data.matrix(y)
    for(j in 1:p){
        ymiss <- is.na(y[,j])
        if(sum(ymiss)>0 & sum(is.na(X[[j]]))>0){
            for(i in 1:kn[j]){
                xmiss<-is.na(X[[j]][, i])
                y[xmiss,j]<-NA
            }
            if(!identical(ymiss,is.na(y[,j])))
                warning("Missing values in X, corresponding elements in y set to NA.")            
        }
    }        
    attributes(y)<-tspy
    kncs<-c(0,cumsum(kn))
    Z<-array(0, dim=c(p, m, n))  
    states<-NULL
    for (j in 1:p){                 
        Z[j,(kncs[j]+1):kncs[j+1] , ] <- t(X[[j]])
        states<-c(states,paste0("beta",1:kn[j],".",j))
    }
    
    #H
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
    T<-diag(m)
    dim(T)<-c(m,m,1)
    
    if(is.null(Q)){
        Q <- array(0, dim=c(m, m, 1))  
        r<-as.integer(m)        
    } else {   
        r<-as.integer(max(dim(Q)[1],1))        
        Q <- array(Q, dim=c(r, r,  (n-1) * (max(dim(Q)[3], 0,na.rm=TRUE) > 1) + 1))        
    } 
    
    R <- array(0, dim=c(m, r, 1))
    R[,,1]<-diag(m)[,1:r,drop=FALSE]         
    
    
    if(!identical(distribution,"Gaussian")){
        if(is.null(u))
            u<-rep(1,n)    
        u<-array(u,dim=n)    
    }
    a1<-matrix(0,m,1)
    rownames(a1)<-states
    object<-list(y=y,Z=Z,H=H,T=T,R=R,Q=Q,a1=a1,P1=matrix(0,m,m),P1inf=diag(m),u=as.double(u), p=p,n=n,m=m,k=r,distribution=distribution,H_type=H_type,tolF=tolF,tol0=tol0)
    class(object) <- c("SSModel","regSSM")
    if (transform %in% c("ldl", "augment") & sum(is.na(H))==0) 
        object <- transformSSM(object, type = transform)
    
    invisible(object)
} 
