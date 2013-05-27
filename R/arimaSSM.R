#' Create a State Space Representation of ARIMA Model
#'
#' Function \code{arimaSSM} creates a state space representation of ARIMA model.
#' 
#' The linear Gaussian state space model is given by
#' 
#' \deqn{y_t = Z_t \alpha_t + \epsilon_t,}{y[t] = Z[t]\alpha[t] + \epsilon[t], (observation equation)}
#' 
#' \deqn{\alpha_{t+1} = T_t \alpha_t + R_t \eta_t,}{\alpha[t+1] = T[t]\alpha[t] + R[t]\eta[t], (transition equation)}
#' 
#' where \eqn{\epsilon_t ~ N(0,H_t)}{\epsilon[t] ~ N(0,H[t])}, \eqn{\eta_t ~ N(0,Q_t)}{\eta[t] ~ N(0,Q[t])} 
#' and \eqn{\alpha_1 ~ N(a_1,P_1)}{\alpha[1] ~ N(a[1],P[1])} independently of each other. In case of non-Gaussian observations,
#' the observation equation is of form \eqn{p(y_t|\theta_t) = p(y_t|Z_t\alpha_t)}{p(y[t]|\theta[t]) = p(y[t]|Z[t]\alpha[t])},
#' with \eqn{p(y_t|\theta_t)}{p(y[t]|\theta[t])} being one of the following:
#'
#' If observations are Poisson distributed, parameter of Poisson distribution
#' is \eqn{u_t\lambda_t}{u[t]\lambda[t]} and \eqn{\theta_t = log(\lambda_t)}{\theta[t]=log(\lambda[t])}. 
#'
#' If observations are from binomial distribution, \eqn{u} is a vector
#' specifying number the of trials at times \eqn{1,\ldots,n}, and \eqn{\theta_t =
#' log[\pi_t/(1-\pi_t)]}{\theta[t] = log(\pi[t]/(1-\pi[t]))}, where \eqn{\pi_t}{\pi[t]} is the probability of success at time \eqn{t}.
#'
#' For non-Gaussian models \eqn{u_t=1}{u[t]=1} as a default. For Gaussian models, parameter is omitted.
#'
#' Only univariate observations are supported when observation equation is non-Gaussian.
#'
#' @export
#' @inheritParams SSModel
#' @seealso \code{\link{regSSM}} for state space representation of a regression model, \code{\link{structSSM}} for structural time series model, and \code{\link{SSModel}} for custom \code{SSModel} object.
#' @param arima A list or a list of lists with components \code{ar}, \code{ma} and \code{d}, giving the autoregression and moving average coefficients, and the degree of differencing for each series. 
#'  If arima is a single list, it is assumed that all \eqn{p} series have same ARIMA structure. Otherwise first sublist gives the ARIMA structure of the first series etc.
#' @param H A \eqn{p \times p}{p*p} covariance matrix (or \eqn{p \times p \times n}{p*p*n} array in of time-varying case) of the disturbance terms 
#'  \eqn{\epsilon_t}{\epsilon[t]} of the observation equation. Default gives \eqn{p \times p}{p*p} zero matrix ie. ordinary ARIMA model without additional noise.
#'  Omitted in case of non-Gaussian distributions. Augment the state vector if you to add want additional noise.
#' @param Q A \eqn{p \times p}{p*p} covariance matrix of the disturbance terms 
#'  \eqn{\eta_t}{\eta[t]} of the system equation. Default is \eqn{p \times p}{p*p} identity matrix ie. ordinary ARIMA model with disturbance terms having unit variance.
arimaSSM <- function(y, arima, H = NULL, Q = NULL, u = NULL, 
        distribution = c("Gaussian", "Poisson", "Binomial"), 
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
    if(!(sum(sapply(arima,is.list))==p)){               
        
        arima2<-vector("list",p)
        for(i in 1:p)
            arima2[[i]]<-arima
        arima<-arima2
    } 
    
    qn<-pn<-dn<-mp<-numeric(p)
    
    for(i in 1:p){
        pn[i] <- length(arima[[i]]$ar) #number of autoregressive parameters
        qn[i] <- length(arima[[i]]$ma) #number of moving average parameters 
        dn[i] <- max(arima[[i]]$d,0)
        mp[i]<-max(pn[i], qn[i] + 1) + dn[i]
    }
    
    
    m<-as.integer(sum(mp))
    
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
        Q <- diag(p)
        dim(Q)<-c(p,p,1)
    } else {            
        Q <- array(Q, dim=c(p, p, 1))        
    } 
    r<-as.integer(p)
    
    
    R <- array(0, dim=c(m, r, 1))
    T <- array(0, dim=c(m, m, 1))
    Z<-array(0, dim=c(p, m, 1))
    P1inf<-P1<-matrix(0, m, m)
    mpcs<-c(0,cumsum(mp))
    for(i in 1:p){
        
        
        Z1<-matrix(0,1,mp[i])
        T1<-matrix(0,mp[i],mp[i])    
        R1<-rep(0,mp[i])
        
        Z1[1,1:(dn[i]+1)]<-1         
        if(dn[i]>0)
        {      
            T1[1:dn[i],1:dn[i]][upper.tri(T1[1:dn[i],1:dn[i]],diag=TRUE)]<-1
            T1[1:dn[i],(dn[i]+1)]<-1
        }   
        if(pn[i]>0)
            T1[(dn[i]+1):(dn[i]+pn[i]),dn[i]+1]<-arima[[i]]$ar
        if(mp[i]>(dn[i]+1))
            T1[(dn[i]+1):(mp[i]-1),(dn[i]+2):mp[i]]<-diag(1,max(pn[i],qn[i]+1)-1)        
        
        R1[dn[i]+1]<-1
        if(qn[i]>0) R1[(dn[i]+2):(dn[i]+1+qn[i])]<-arima[[i]]$ma
        
        Z[i,(mpcs[i]+1):mpcs[i+1],1]<-Z1 
        T[(mpcs[i]+1):mpcs[i+1],(mpcs[i]+1):mpcs[i+1],1]<-T1
        R[(mpcs[i]+1):mpcs[i+1],i,1]<-R1 
        if(dn[i]>0)
            P1inf[(mpcs[i]+1):(mpcs[i]+dn[i]),(mpcs[i]+1):(mpcs[i]+dn[i])]<-diag(1,dn[i])
    }
    
    
    if(sum(dn)>0){
        diffuse<-which(diag(P1inf)==1)        
        P1[-diffuse,-diffuse] <- solve(a = diag((m-length(diffuse))^2) - matrix(kronecker(T[-diffuse,-diffuse,1], T[-diffuse,-diffuse,1]),(m-length(diffuse))^2,(m-length(diffuse))^2), 
                                    b = c(matrix(R[-diffuse,,1],m-length(diffuse),r) %*%Q[,,1]%*%t(matrix(R[-diffuse,,1],m-length(diffuse),r))))
    } else {
        P1[] <- solve(a = diag(m^2) - matrix(kronecker(T[,,1], T[,,1]),m^2,m^2), b = c(matrix(R[,,1],m,r) %*%Q[,,1]%*%t(matrix(R[,,1],m,r))))       
    }
    
    if(!identical(distribution,"Gaussian")){
        if(is.null(u))
            u<-rep(1,n)    
        u<-array(u,dim=n)    
    }
    a1<-matrix(0,m,1)
    states<-NULL
    for (j in 1:p){                         
        states<-c(states,paste0("state",1:mp[j],".",j))
    }
    
    rownames(a1)<-states
    object<-list(y=y,Z=Z,H=H,T=T,R=R,Q=Q,a1=a1,P1=P1,P1inf=P1inf,u=as.double(u),p=p,n=n,m=m,k=r,distribution=distribution,H_type=H_type,tolF=tolF,tol0=tol0)
    class(object) <- c("SSModel","arimaSSM")
    if (transform %in% c("ldl", "augment") & sum(is.na(H))==0) 
        object <- transformSSM(object, type = transform)
    
    invisible(object)
} 
