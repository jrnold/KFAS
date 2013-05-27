#' Create a Structural Time Series State Space Model
#'
#' Function \code{structSSM} creates a state space representation of structural time series. 
#' 
#' The structural time series model has form
#' 
#' \deqn{y_t = \mu_t + \gamma_t + \epsilon_t, \quad \epsilon_t \sim N(0, H_t) }{y[t] = \mu[t] + \epsilon[t], \epsilon[t] ~ N(0, H[t])}
#' 
#' \deqn{\mu_{t+1} = \mu_t + \nu_t + \xi_t, \quad \xi_t \sim N(0, Q_{level,t})}{\mu[t+1] = \mu[t] + \nu[t] + \xi[t], \xi[t] ~ N(0, Q[level,t])}
#' 
#' \deqn{\nu_{t+1} = \nu_t + \zeta_t, \quad \zeta_t \sim  N(0, Q_{slope,t})}{\nu[t+1] = \nu[t] + \zeta[t], \zeta[t] ~ N(0, Q[slope,t]) }
#' 
#' with seasonal component being either time domain form
#' \deqn{\gamma_{t+1} = -\sum_{j=1}^{s-1}\gamma_{t+1-j} + \omega_t, \quad \omega_t \sim N(0,Q_{seasonal,t}), }{\gamma[t+1] = -\gamma[t] - \ldots - \gamma[t-s+2] + \omega[t], \omega[t] ~ N(0,Q[seasonal,t]), }
#' 
#' or frequency domain form where
#' 
#' \deqn{\gamma_{t} = \sum_{j=1}^{\lfloor s/2 \rfloor}\gamma_{j,t} }{\gamma[t] = \gamma[1,t] + \ldots + \gamma[[s/2],t], }
#' \deqn{\gamma_{j,t+1} = \gamma_{j,t} cos\lambda_j + \gamma^{\ast}_{j,t} sin\lambda_j + \omega_{j,t},  }{\gamma_{j,t+1} = \gamma[j,t] cos\lambda[j] + \gamma*[j,t] sin\lambda[j] + \omega[j,t],  }
#' \deqn{\gamma^{\ast}_{j,t+1} = - \gamma_{j,t} sin\lambda_j + \gamma^{\ast}_{j,t} cos\lambda_j + \omega^\ast_{j,t}, j=1,\ldots, \lfloor s/2 \rfloor,  }{\gamma*[j,t+1] = - \gamma[j,t] sin\lambda[j] + \gamma*[j,t] cos\lambda[j] + \omega*[j,t], j=1,..., [s/2],  }
#' with \eqn{\omega_{j,t}}{\omega[j,t]} and \eqn{\omega^\ast_{j,t}}{\omega*[j,t]} being independently distributed variables with \eqn{N(0, Q_{seasonal,t})}{N(0, Q[seasonal,t])} distribution and \eqn{\lambda_j = 2\pi j/s}{\lambda[j] = 2\pi j/s}.
#' 
#' Explanatory variables can also be added to the model; in \code{structSSM} function it is assumed that same explanatory variables are used for all series.
#' See \code{\link{regSSM}} and \code{\link{+}} for more complicated settings.
#' @export
#' @inheritParams SSModel
#' @seealso \code{\link{arimaSSM}} for state space representation of ARIMA model, \code{\link{regSSM}} for state space representation of a regression model, \code{\link{SSModel}} for custom \code{SSModel} object and \code{\link{KFAS}} for general information regarding the package and examples of its usage.
#' @param trend A character vector defining the type of the level component of the model. 
#' For multivariate series, either one type, it is assumed that all \eqn{p} series have same type of trend components. 
#' Possible values are \code{"level"} (local level model) and \code{"slope"} (local linear trend model). Default is \code{"level"}.
#' @param seasonal A character vector defining the type of the seasonal component of the model. 
#' For multivariate series, it is assumed that all \eqn{p} series have same type of seasonal components.
#'  Possible values are \code{"none"} (no seasonal), \code{"time"} (time domain form) and \code{"time"} (frequency domain form).
#' The length of the seasonal pattern is taken as the frequency attribute of the time series object \code{y}. Default is \code{"none"}.
#' @param X A \eqn{n \times k}{n*k} matrix of explanatory variables, with each column containing one explanatory variable. It is assumed that all \eqn{p} 
#' series use same explanatory variables.  
#' @param H A \eqn{p \times p}{p*p} covariance matrix (or \eqn{p \times p \times n}{p*p*n} array in time-varying case) of the disturbance terms 
#'   \eqn{\epsilon_t}{\epsilon[t]} of the observation equation. Default gives \eqn{p \times p}{p*p} zero matrix.
#'   Omitted in case of non-Gaussian distributions. Augment the state vector if you want to add additional noise.
#' @param Q.level A scalar or \eqn{p \times p}{p*p} covariance matrix (or \eqn{p \times p \times n}{p*p*n} array in time-varying case) of the disturbance terms \eqn{\xi_t}{\xi[t]} corresponding to the level process \eqn{\mu_t}{\mu[t]}. 
#'  Default gives diagonal matrix with NA's on diagonal. 
#' @param Q.slope A scalar or \eqn{p \times p}{p*p} covariance matrix (or \eqn{p \times p \times n}{p*p*n} array in time-varying case) of the disturbance terms \eqn{\zeta_t}{\xi[t]} corresponding to the slope process \eqn{\nu_t}{\nu[t]}. 
#'  Default gives diagonal matrix with NA's on diagonal. Omitted if \code{trend="level"}. 
#' @param Q.seasonal scalar or A \eqn{p \times p}{p*p} covariance matrix (or \eqn{p \times p \times n}{p*p*n} array in time-varying case) of the disturbance terms \eqn{\omega_t}{\omega[t]} corresponding to the seasonal process \eqn{\gamma_t}{\gamma[t]}. 
#'  Default gives diagonal matrix with NA's on diagonal. Omitted if \code{seasonal="none"}.  There are several \eqn{\omega_t}{\omega[t]} processes in the frequency domain case, but they are all identically distributed, so only the (co)variance structure of one of them need to be defined.
#' @param Q.regression A scalar or \eqn{xn \times xn}{xn*xn} covariance matrix (or \eqn{xn \times xn \times n}{xn*xn*n} array in time-varying case) of 
#' the disturbance terms corresponding to the regression coefficient processes. 
#'  Default gives zero matrix i.e. ordinary time-invariant regression. 
structSSM <- function(y, trend="level", seasonal="none", X=NULL, H=NULL, Q.level=NULL, Q.slope=NULL, Q.seasonal=NULL, Q.regression=NULL, u=NULL,
        distribution = c("Gaussian", "Poisson", "Binomial"), 
        transform = c("none", "ldl", "augment"), tolF = .Machine$double.eps^0.5, tol0 = .Machine$double.eps^0.5) {
    
    transform <- match.arg(arg=transform, choices = c("none", "ldl", "augment"))    
    distribution <- match.arg(arg=distribution, choices = c("Gaussian", "Poisson", "Binomial"))    
    trend <- match.arg(arg=trend, choices = c("level", "slope")) 
    seasonal <- match.arg(arg=seasonal, choices = c("none","time", "frequency")) 
    
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
    
    
    s<-frequency(y)
    
    tr<-identical(trend,"slope")
    
    m <- k<- p+p*tr
    
    if(p==1){      
        states<-"level"
        if(tr)
            states<-c(states,paste("slope",sep=""))
    } else{      
        states<-paste("level.",1:p,sep="")
        if(tr)
            states<-c(states,paste("slope.",1:p,sep=""))
    }
    
    
    
    if(!identical(seasonal,"none")){          
        m <- m + p*(s - 1)
        if(identical(seasonal,"time")){
            k <- k + p  
            if(p==1){
                states<-c(states,paste(c(paste("seasonal"),
                                        paste("seasonal_",-1:(-s+2),sep="")),sep=""))
            } else{
                states<-c(states,paste(c(rep(paste("seasonal."),each=p),
                                        rep(paste("seasonal_",-1:(-s+2),".",sep=""),each=p)),1:p,sep=""))
            }
            
        }else {                   
            k <- k + p*(s - 1)
            if(p==1){
                states<-c(states,paste(rep(c("seasonal_","seasonal*_"),length.out=(s - 1)),rep(1:floor(s/2),each=2,length.out=(s - 1)),sep=""))
            } else{
                states<-c(states,paste(rep(c("seasonal_","seasonal*_"),each=p,length.out=p*(s - 1)),rep(1:floor(s/2),each=2*p,length.out=p*(s - 1)),".",1:p,sep=""))
            }
            
        }
    }
    
    
    
    
    
    
    
    
    
    
    Z <- matrix(0, p, m)
    T <- matrix(0, m, m)
    R <- matrix(0, m, k)
    
    diag(Z)[1:p]<-diag(T)[1:p]<-diag(R)[1:p]<-1
    
    if(tr){
        T[(p+1):(2*p),(p+1):(2*p)][1 + 0:(p - 1) * (p + 1)]<-1
        R[(p+1):(2*p),(p+1):(2*p)][1 + 0:(p - 1) * (p + 1)]<-1       
        T[1:p,(p+1):(2*p)][1 + 0:(p - 1) * (p + 1)]<-1
        
    }
    
    if(!identical(seasonal,"none")){    
        if(identical(seasonal,"time")){
            if(s<2) stop("Length of the seasonal must be at least 3.")
            Z[1:p,(p+p*tr+1):(2*p+p*tr)][1 + 0:(p - 1) * (p + 1)]<-1
            for(i in 1:p)
                T[p+p*tr+i,(p+p*tr+i-1+seq(1,p*(s-1),by=p))] <- -1
            #T[(p+p*tr+1):(2*p+p*tr),(p+p*tr+1):m][row(matrix(0,p,p*(s-1))) <= col(matrix(0,p,p*(s-1)))] <- -1
            diag(T[(2*p+p*tr+1):m,(p+p*tr+1):m])<- 1            
            R[(p+p*tr+1):k,(p+p*tr+1):k][1 + 0:(k-(p+p*tr) - 1) * (k-(p+p*tr) + 1)]<-1
        } else {
            if(s<3) stop("Length of the seasonal must be at least 3.")
            
            Z[1:p,(p+p*tr+1):m]<-rep(c(diag(p),rep(0,p^2)),length.out=p*p*(s-1))
            
            for(j in 1:floor((s-1)/2)){
                lambda<-2*pi*j/s
                T[(p+p*tr+2*(j-1)*p+1):(p+p*tr+2*j*p),(p+p*tr+2*(j-1)*p+1):(p+p*tr+2*j*p)][1 + 0:(2*p - 1) * (2*p + 1)]<-cos(lambda)
                T[(p+p*tr+2*(j-1)*p+1):(p+p*tr+2*(j-1)*p+p),p+(p+p*tr+2*(j-1)*p+1):(p+p*tr+2*(j-1)*p+p)][1 + 0:(p - 1) * (p + 1)]<-sin(lambda)
                T[p+(p+p*tr+2*(j-1)*p+1):(p+p*tr+2*(j-1)*p+p),(p+p*tr+2*(j-1)*p+1):(p+p*tr+2*(j-1)*p+p)][1 + 0:(p - 1) * (p + 1)]<- - sin(lambda)                  
            }
            if(s%%2==0)
                diag(T)[(m-p+1):m]<- -1
            R<-diag(k)
            
            
        }
    }   
    
    dim(Z)<-c(p,m,1)
    dim(T)<-c(m,m,1)
    dim(R)<-c(m,k,1)
    
    
    if(!is.null(X)){                   
        if(!is.ts(X) & (!identical(dim(X)[1],n) | is.null(dim(X)[2]))) stop("X is not an n times k matrix or a time series object.")
        X<-data.matrix(X)
        kn <- dim(X)[2]
        
        tspy<-attributes(y)
        y<-data.matrix(y)
        ymiss <- is.na(y)
        if(sum(ymiss)>0 & sum(is.na(X))>0){
            xmiss<-is.na(rowSums(X))
            y[xmiss,]<-NA
            if(!identical(ymiss,is.na(y)))
                warning("Missing values in X, corresponding elements in y set to NA.")            
        }
        attributes(y)<-tspy       
        m2 <- m + kn*p
        Z2<-array(0, dim=c(p, m2, n))   
        Z2[,1:m,]<-Z
        for(i in 1:kn)
            for(j in 1:p)
                Z2[j,(m+j+p*(i-1)), ] <- X[,i]
        Z<-Z2
        R2<-array(0, dim=c(m2, k+kn*p, 1))
        R2[1:m,1:k,1]<-R
        R2[(m+1):m2,(k+1):(k+kn*p),1][1 + 0:(kn*p - 1) * (kn*p + 1)]<-1
        R<-R2
        T2<-array(0, dim=c(m2, m2, 1))
        T2[1:m,1:m,1]<-T
        T2[(m+1):m2,(m+1):m2,1][1 + 0:(kn*p - 1) * (kn*p + 1)]<-1
        T<-T2
        
        m<-m2
        k<-k+kn*p
        
        
        states<-c(states,paste(rep(paste("beta",1:kn,".",sep=""),each=p),1:p,sep=""))
    } else kn<-0
    
    storage.mode(m)<- storage.mode(k)<-"integer"
    
    
    if(distribution!="Gaussian"){
        H_type<-"Omitted"
        H <- array(0, dim=c(p, p, 1))
    } else {
        if(is.null(H)){
            H <- diag(NA,p)
            dim(H)<-c(p,p,1)
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
    
    tv<-max(dim(Q.level)[3], (trend=="slope")*dim(Q.slope)[3],(seasonal!="none")*dim(Q.seasonal)[3],(seasonal!="none")*dim(Q.regression)[3],0,na.rm=TRUE)
    Q <- array(0, dim=c(k, k,  (n-1) * (tv > 1) + 1))      
    if(!is.null(Q.level)){
        #if(dim(Q.level)[1:2]!=c(k,k))
        Q[1:p,1:p,]<-Q.level               
    } else {
        Q[1:p,1:p,]<-diag(NA,p)
    }
    if(tr){
        if(!is.null(Q.slope)){      
            Q[(p+1):(2*p),(p+1):(2*p),]<-Q.slope               
        } else {
            Q[(p+1):(2*p),(p+1):(2*p),]<-diag(NA,p)
        }
    }
    if(seasonal!="none"){
        
        if(!is.null(Q.seasonal)){    
            if(seasonal=="time"){
                Q[(p+tr*p+1):(2*p+tr*p),(p+tr*p+1):(2*p+tr*p),]<-Q.seasonal    
            } else {
                for(j in 1:(s-1))
                    Q[(p+tr*p+1+p*(j-1)):(p+tr*p+j*p),(p+tr*p+1+p*(j-1)):(p+tr*p+j*p),]<-Q.seasonal
            }
        } else {
            Q[(p+tr*p+1):(k-p*kn),(p+tr*p+1):(k-p*kn),][1 + 0:(k-p*kn-p+tr*p - 1) * (k-p*kn-p+tr*p + 1)]<-NA
        } 
    } 
    if(!is.null(X) & !is.null(Q.regression)){
        for(i in 1:p)
            Q[(k-kn*p+1+(i-1)*p):(k-kn*p+(i-1)*p+kn),(k-kn*p+1+(i-1)*p):(k-kn*p+(i-1)*p+kn),]<-Q.regression   
    }
    
    
    
    if(!identical(distribution,"Gaussian")){
        if(is.null(u))
            u<-rep(1,n)    
        u<-array(u,dim=n)    
    }
    a1<-matrix(0,m,1)
    rownames(a1)<-states
    object<-list(y=y,Z=Z,H=H,T=T,R=R,Q=Q,a1=a1,P1=matrix(0,m,m),P1inf=diag(m),u=as.double(u),p=p,n=n,m=m,k=k,distribution=distribution,H_type=H_type,tolF=tolF,tol0=tol0)
    class(object) <- c("SSModel","structSSM")
    if (transform %in% c("ldl", "augment") & sum(is.na(H))==0) 
        object <- transformSSM(object, type = transform)
    
    invisible(object)
    
} 
