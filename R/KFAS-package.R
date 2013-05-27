#' KFAS: Functions for Gaussian and Non-Gaussian State Space Models
#'
#' Package KFAS contains functions for Kalman filtering, smoothing and simulation of linear state space models
#' with exact diffuse initialization.
#'
#' The linear Gaussian state space model is given by
#' 
#' \deqn{y_t = Z_t \alpha_t + \epsilon_t,}{y[t] = Z[t]\alpha[t] + \epsilon[t], (observation equation)}
#' 
#' \deqn{\alpha_{t+1} = T_t \alpha_t + R_t \eta_t,}{\alpha[t+1] = T[t]\alpha[t] + R[t]\eta[t], (transition equation)}
#' 
#' where \eqn{\epsilon_t \sim N(0,H_t)}{\epsilon[t] ~ N(0,H[t])}, \eqn{\eta_t \sim N(0,Q_t)}{\eta[t] ~ N(0,Q[t])} 
#' and \eqn{\alpha_1 \sim N(a_1,P_1)}{\alpha[1] ~ N(a[1],P[1])} independently of each other.
#' 
#' All system and covariance matrices Z, H, T, R and Q can be time-varying, and partially or totally missing observations \eqn{y_t}{y[t]} are allowed.
#' 
#' Covariance matrices H and Q has to be positive semidefinite. 
#' 
#' Dimensions of system matrices are
#' 
#' \tabular{rl}{
#'  Z \tab \eqn{p \times m \times 1}{p*m*1} or \eqn{p \times m \times n}{p*m*n} in time varying case \cr
#'  H \tab \eqn{p \times p \times 1}{p*p*1} or \eqn{p \times p \times n}{p*p*n} in time varying case \cr
#'  T \tab \eqn{m \times m \times 1}{m*m*1} or \eqn{m \times m \times n}{m*m*n} in time varying case \cr
#'  R \tab \eqn{m \times k \times 1}{m*k*1} or \eqn{m \times k \times n}{m*k*n} in time varying case \cr
#'  Q \tab \eqn{k \times k \times 1}{k*k*1} or \eqn{k \times k \times n}{k*k*n} in time varying case \cr
#'  }
#' 
#' In case of non-Gaussian observations, the observation equation is of form 
#' \deqn{p(y_t|\theta_t) = p(y_t|Z_t\alpha_t)}{p(y[t]|\theta[t]) = p(y[t]|Z[t]\alpha[t]),}
#' with \eqn{p(y_t|\theta_t)}{p(y[t]|\theta[t])} being one of the following:
#'
#' If observations \eqn{y_t}{y[t]} are Poisson distributed, parameter of Poisson distribution
#' is \eqn{u_t\lambda_t}{u[t]\lambda[t]} and \eqn{\theta_t = log(\lambda_t)}{\theta[t]=log(\lambda[t])}, 
#' where \eqn{u_t}{u[t]} is so called offset term.
#' 
#' 
#' If observations eqn{y_t}{y[t]} are from binomial distribution, \eqn{u} is a vector
#' specifying number the of trials at times \eqn{1,\ldots,n}, and \eqn{\theta_t =
#' log[\pi_t/(1-\pi_t)]}{\theta[t] = log(\pi[t]/(1-\pi[t]))}, where \eqn{\pi_t}{\pi[t]} is the probability of success at time \eqn{t}.
#'
#' For non-Gaussian models \eqn{u_t=1}{u[t]=1} as a default. For Gaussian models, parameter is omitted.
#'
#' Only univariate observations are supported when observation equation is non-Gaussian.
#'
#' For the unknown elements of initial state vector \eqn{a_1}{a[1]}, KFS
#' uses exact diffuse initialization by Koopman and Durbin (2000, 2001, 2003), where the unknown initial states are set to have a zero mean and infinite variance, so
#' \deqn{P_1 = P_{\ast,1} + \kappa P_{\infty,1},}{P[1] = P[*,1] + \kappaP[inf,1],}
#' with \eqn{\kappa} going to infinity and \eqn{P_{\infty,1}}{P[inf,1]} being diagonal matrix with ones on diagonal elements corresponding to unknown initial states. 
#' 
#' Diffuse phase is continued until rank of \eqn{P_{\infty,t}}{P[inf,t]} becomes zero. Rank of \eqn{P_{\infty}}{P[inf]} decreases by 1, if \eqn{F_\infty>tolF>0}{F[inf]>tolF>0}. 
#' Usually the number of diffuse time points equals the number unknown elements of initial state vector, but missing observations or time-varying Z can affect this.
#' See Koopman and Durbin (2000, 2001, 2003) for details for exact diffuse and non-diffuse filtering.
#'
#' To lessen the notation and storage space, KFAS uses letters P, F and K for non-diffuse part of the corresponding matrices, omitting the asterisk in diffuse phase.
#' 
#' All functions of KFAS use the univariate approach (also known as sequential processing, see Anderson and Moore (1979))
#' which is from Koopman and Durbin (2000, 2001). In univariate approach the observations are introduced one element at the time.
#' Therefore the prediction error variance matrices F and Finf does not need to be non-singular, as
#' there is no matrix inversions in univariate approach algorithm.  This provides more stable and possibly more faster filtering and smoothing than normal multivariate Kalman
#' filter algorithm. If covariance matrix H is not diagonal, it is possible to transform the model by either using LDL decomposition on H, or augmenting the state vector with \eqn{\epsilon} disturbances. See \code{\link{transformSSM}} for more details.
#'
#'
#' @references Koopman, S.J. and Durbin J. (2000).  Fast filtering and
#' smoothing for non-stationary time series models, Journal of American
#' Statistical Assosiation, 92, 1630-38.
#'
#' Koopman, S.J. and Durbin J. (2001).  Time Series Analysis by State Space
#' Methods. Oxford: Oxford University Press.
#'
#' Koopman, S.J. and Durbin J. (2003).  Filtering and smoothing of state vector
#' for diffuse state space models, Journal of Time Series Analysis, Vol. 24,
#' No. 1.
#' #' Shumway, Robert H. and Stoffer, David S. (2006).  Time Series Analysis and
#' Its Applications: With R examples.  \cr
#' @docType package
#' @name KFAS
#' @aliases KFAS
#' @examples
#'
#' library(KFAS)
#'
#' # Example of local level model for Nile series
#'
#' y<-Nile
#' modelNile<-structSSM(y=y)
#'
#' fit<-fitSSM(inits=c(0.5*log(var(Nile)),0.5*log(var(Nile))),model=modelNile)
#' # Filtering and state smoothing
#' kfsNile<-KFS(fit$model,smoothing="state") 
#' # Simple plot of series and the smoothed signal = Z*alphahat
#' plot(kfsNile,col=1:2)
#'
#'
#' # Confidence intervals for the state
#' lows<-c(kfsNile$alphahat-qnorm(0.95)*sqrt(c(kfsNile$V)))
#' ups<-c(kfsNile$alphahat+qnorm(0.95)*sqrt(c(kfsNile$V)))
#' plot.ts(cbind(y,c(kfsNile$alphahat),lows,ups), plot.type="single", col=c(1:2,3,3), 
#'         ylab="Predicted Annual flow", main="River Nile")
#'
#'
#' # Missing observations, using same parameter estimates
#'
#' y<-Nile
#' y[c(21:40,61:80)]<-NA
#' modelNile<-structSSM(y=y,H=fit$model$H,Q.level=fit$model$Q)
#'
#' kfsNile<-KFS(modelNile,smoothing="state")
#'
#' # Filtered state
#' plot.ts(cbind(y,c(NA,kfsNile$a[,-c(1,101)])), plot.type="single", col=c(1:2,3,3), 
#'        ylab="Predicted Annual flow", main="River Nile")
#'
#' #  Smoothed state
#' plot.ts(cbind(y,c(kfsNile$alp)), plot.type="single", col=c(1:2,3,3), 
#'        ylab="Predicted Annual flow", main="River Nile")
#'
#' 
#' # Prediction of missing observations
#' predictNile<-predict(kfsNile)
#' lows<-predictNile$y-qnorm(0.95)*sqrt(c(predictNile$F))
#' ups<-predictNile$y+qnorm(0.95)*sqrt(c(predictNile$F))
#'
#' plot.ts(cbind(y,predictNile$y,lows,ups), plot.type="single", col=c(1:2,4,4), 
#'         ylab="Predicted Annual flow", main="River Nile")
#' 
#'
#' 
#' # Example of multivariate local level model with only one state
#' # Two series of average global temperature deviations for years 1880-1987
#' # See Shumway and Stoffer (2006), p. 327 for details
#'
#' data(GlobalTemp) 
#'
#' modelTemp<-SSModel(y=GlobalTemp, Z = matrix(1,nrow=2), T=1, R=1, H=matrix(NA,2,2), 
#'         Q=NA, a1=0, P1=0, P1inf=1)
#'
#' # Estimating the variance parameters
#' 
#' fit<-fitSSM(inits=c(0.5*log(.1),0.5*log(.1),0.5*log(.1),0),model=modelTemp) 
#'
#' outTemp<-KFS(fit$model,smooth="both")
#'
#'
#' ts.plot(cbind(modelTemp$y,t(outTemp$alphahat)),col=1:3)
#' legend("bottomright",legend=c(colnames(GlobalTemp), "Smoothed signal"), col=1:3, lty=1)
#' 
#' 
#' # Example of multivariate series where first series follows stationary ARMA(1,1) 
#' # process and second stationary AR(1) process.
#'
#' y1<-arima.sim(model=list(ar=0.8,ma=0.3), n=100, sd=0.5)
#'
#' model1<-arimaSSM(y=y1,arima=list(ar=0.8,ma=0.3),Q=0.5^2)
#' 
#' y2<-arima.sim(model=list(ar=-0.5), n=100)
#'
#' model2<-arimaSSM(y=y2,arima=list(ar=-0.5))
#'
#' model<-model1+model2
#' 
#' # Another way:
#' 
#' modelb<-arimaSSM(y=ts.union(y1,y2),arima=list(list(ar=0.8,ma=0.3),list(ar=-0.5)),Q=diag(c(0.5^2,1)))
#'
#' f.out<-KFS(model)
#' 
#' 
#' # Drivers
#' model<-structSSM(y=log(Seatbelts[,"drivers"]),trend="level",seasonal="time",
#'        X=cbind(log(Seatbelts[,"kms"]),log(Seatbelts[,"PetrolPrice"]),Seatbelts[,c("law")]))
#' fit<-fitSSM(inits=rep(-1,3),model=model)
#' out<-KFS(fit$model,smoothing="state")
#'
#'
#' plot(out,lty=1:2,col=1:2,main="Observations and smoothed signal with and without seasonal component")
#' lines(signal(out,states=c(1,13:15))$s,col=4,lty=1)
#' legend("bottomleft",legend=c("Observations", "Smoothed signal","Smoothed level"), col=c(1,2,4), lty=c(1,2,1))
#'
#'
#' # Multivariate model with constant seasonal pattern in frequency domain
#'
#' model<-structSSM(y=log(Seatbelts[,c("front","rear")]),trend="level",seasonal="freq",
#'        X=cbind(log(Seatbelts[,c("kms")]),log(Seatbelts[,"PetrolPrice"]),Seatbelts[,"law"]),
#'        H=NA,Q.level=NA,Q.seasonal=0)
#'
#' sbFit<-fitSSM(inits=rep(-1,6),model=model)
#'
#' out<-KFS(sbFit$model,smoothing="state")
#'
#' ts.plot(signal(out,states=c(1:2,25:30))$s,model$y,col=1:4)
#'
#' # Poisson model
#' model<-structSSM(y=Seatbelts[,"VanKilled"],trend="level",seasonal="time",X=Seatbelts[,"law"],
#'        distribution="Poisson")
#' 
#' # Estimate variance parameters
#' fit<-fitSSM(inits=rep(0.5*log(0.005),2), model=model)  
#'
#' model<-fit$model
#' # Approximating model, gives also the conditional mode of theta
#' amod<-approxSSM(model)
#' out.amod<-KFS(amod,smoothing="state")
#' 
#' # State smoothing via importance sampling
#' out<-KFS(model,nsim=1000)
#'
#' # Observations with exp(smoothed signal) computed by 
#' # importance sampling in KFS, and by approximating model
#' ts.plot(cbind(model$y,out$yhat,exp(amod$theta)),col=1:3) # Almost identical
#'
#' # It is more interesting to look at the smoothed values of exp(level + intervention)
#' lev1<-exp(signal(out,states=c(1,13))$s)
#' lev2<-exp(signal(out.amod,states=c(1,13))$s)
#' # These are slightly biased as E[exp(x)] > exp(E[x]), better to use importance sampling:
#' vansample<-importanceSSM(model,save.model=FALSE,nsim=250) 
#' # nsim is number of independent samples, as default two antithetic variables are used, 
#' # so total number of samples is 1000.
#'
#' w<-vansample$weights/sum(vansample$weights)
#' level<-colSums(t(exp(vansample$states[1,,]+model$Z[1,13,]*vansample$states[13,,]))*w)
#' ts.plot(cbind(model$y,lev1,lev2,level),col=1:4) #' Almost identical results
#'
#' # Confidence intervals (no seasonal component)
#' 
#' varlevel<-colSums(t(exp(vansample$states[1,,]+model$Z[1,13,]*vansample$states[13,,])^2)*w)-level^2 
#' intv<-level + qnorm(0.975)*varlevel%o%c(-1,1)
#' ts.plot(cbind(model$y,level,intv),col=c(1,2,3,3))
#' 
#' # Simulation error
#' 
#' # Mean estimation error of the single draw with 2 antithetics
#' level2<-t(exp(vansample$states[1,,1:250]+model$Z[1,13,]*vansample$states[13,,1:250])-level)*w[1:250]+
#'   t(exp(vansample$states[1,,251:500]+model$Z[1,13,]*vansample$states[13,,251:500])-level)*w[251:500]+
#'   t(exp(vansample$states[1,,501:750]+model$Z[1,13,]*vansample$states[13,,501:750])-level)*w[501:750]+
#'   t(exp(vansample$states[1,,751:1000]+model$Z[1,13,]*vansample$states[13,,751:1000])-level)*w[751:1000]
#'      
#' varsim<-colSums(level2^2)
#' ts.plot(sqrt(varsim/varlevel)*100)
#' 
#' # Same without antithetic variables
#' vansamplenat<-importanceSSM(model,save.model=FALSE,nsim=1000,antithetics=FALSE)
#' w<-vansamplenat$weights/sum(vansamplenat$weights)
#' levelnat<-colSums(t(exp(vansamplenat$states[1,,]+
#'           model$Z[1,13,]*vansamplenat$states[13,,]))*w)     
#' varsimnat<-colSums((t(exp(vansamplenat$states[1,,]+
#'             model$Z[1,13,]*vansamplenat$states[13,,])-levelnat)*w)^2)
#' varlevelnat<-colSums(t(exp(vansamplenat$states[1,,]+
#'             model$Z[1,13,]*vansamplenat$states[13,,])^2)*w)-levelnat^2 
#' ts.plot(sqrt(varsimnat/varlevelnat)*100) 
#' ts.plot((sqrt(varsimnat)-sqrt(varsim))/sqrt(varsimnat)*100) 

NULL
#' Oxford-Cambridge boat race results 1829-2000
#'
#' Results of the annual boat race between universities of Oxford (0) and Cambridge (1).
#'
#' @name boat
#' @docType data
#' @format An time series object containing 172 observations.
#' @references  Koopman, S.J. and Durbin J. (2001).  Time Series Analysis by State Space Methods. Oxford: Oxford University Press.
#' @source http://www.ssfpack.com/DKbook.html
#' @keywords datasets

NULL
#' Two series of average global temperature deviations for years 1880-1987
#'
#' This data set contains two series of average global temperature deviations
#' for years 1880-1987. These series are same as used in Shumway and Stoffer
#' (2006), where they are known as HL and Folland series. For more details, see
#' Shumway and Stoffer (2006, p. 327).
#'
#'
#' @name GlobalTemp
#' @docType data
#' @format An time series object containing 108 times 2 observations.
#' @references Shumway, Robert H. and Stoffer, David S. (2006). Time Series
#' Analysis and Its Applications: With R examples.
#' @source http://lib.stat.cmu.edu/general/stoffer/tsa2/
#' @keywords datasets
NULL
#' Combine State Space Model Objects of class \code{SSModel}
#' 
#' Second model needs to have either only duplicate time series with first model, or no identical series at all.
#' @S3method + SSModel
#' @method + SSModel
#' @param e1 ,
#' @param e2 Models to be combined.
#' @return \item{model}{Combined model.}
#' @rdname combine.SSModel
"+.SSModel" <- function(e1, e2) {    
    
    if(!identical(e1$n,e2$n))
        stop("Models cannot be combined as the number of time points differ.")
    if(!identical(e1$distribution,e2$distribution))      
        stop("Models cannot be combined as their distributions differ.")
    
    
    if(!identical(e1$H_type,e2$H_type)){
        if("Diagonal"%in%c(e1$H_type,e2$H_type)){
            if(!("Augmented"%in%c(e1$H_type,e2$H_type))){
                H_type<-c(e1$H_type,e2$H_type)[-which("Diagonal"==c(e1$H_type,e2$H_type))]    
            } else H_type<-"Diagonal"            
        } else {
            stop("Model with augmented H can only be combined with model with augmented or diagonal H.")
        }        
    } else {
        H_type <- e1$H_type        
    }  
    
    distribution <- e1$distribution
    
    y<-ts.union(e1$y,e2$y)
    
    dups<-duplicated(t(y))    
    sameSeries<-sum(dups)==e1$p & e1$p==e2$p
    differentSeries<-sum(dups)==0
    
    if(sameSeries | differentSeries){         
        p<-e1$p+differentSeries*e2$p        
        m<-as.integer(e1$m+e2$m)
        k<-as.integer(e1$k+e2$k)
        n<-e1$n
        Z<-array(0, dim=c(p, m, (n-1) * (max(dim(e1$Z)[3], dim(e2$Z)[3], 0) > 1) + 1))
        H<-array(0, dim=c(p, p, (n-1) * (max(dim(e1$H)[3], dim(e2$H)[3], 0) > 1) + 1))
        T<-array(0, dim=c(m, m, (n-1) * (max(dim(e1$T)[3], dim(e2$T)[3], 0) > 1) + 1))
        R<-array(0, dim=c(m, k, (n-1) * (max(dim(e1$R)[3], dim(e2$R)[3], 0) > 1) + 1))
        Q<-array(0, dim=c(k, k, (n-1) * (max(dim(e1$Q)[3], dim(e2$Q)[3], 0) > 1) + 1))
        
        Z[1:e1$p,1:e1$m,]<-e1$Z
        H[1:e1$p,1:e1$p,]<-e1$H
        T[1:e1$m,1:e1$m,]<-e1$T
        R[1:e1$m,1:e1$k,]<-e1$R
        Q[1:e1$k,1:e1$k,]<-e1$Q
        P1<-P1inf<-matrix(0,m,m)
        P1[1:e1$m,1:e1$m]<-e1$P1
        P1[(e1$m+1):m,(e1$m+1):m]<-e2$P1
        P1inf[1:e1$m,1:e1$m]<-e1$P1inf
        P1inf[(e1$m+1):m,(e1$m+1):m]<-e2$P1inf
        a1<-matrix(0,m,1)
        a1[1:e1$m,1]<-e1$a1
        a1[(e1$m+1):m,1]<-e2$a1
        rownames(a1)<-c(rownames(e1$a1),rownames(e2$a1))
        tol0<-max(e1$tol0,e2$tol0)
        tolF<-max(e1$tolF,e2$tolF)
        
        if(sameSeries){  
            y<-e1$y
            Z[1:p,(e1$m+1):m,]<-e2$Z
            H <- H + e2$H              
            T[(e1$m+1):m,(e1$m+1):m,]<-e2$T
            R[(e1$m+1):m,(e1$k+1):k,]<-e2$R
            Q[(e1$k+1):k,(e1$k+1):k,]<-e2$Q
            
        } else {
            Z[(e1$p+1):p,(e1$m+1):m,]<-e2$Z
            H[(e1$p+1):p,(e1$p+1):p,]<-e2$H
            T[(e1$m+1):m,(e1$m+1):m,]<-e2$T
            R[(e1$m+1):m,(e1$k+1):k,]<-e2$R
            Q[(e1$k+1):k,(e1$k+1):k,]<-e2$Q
            
        }
    } else {
        stop("Only models with either identical or completely non-identical time series y can be combined. ")        
    }   
    
    
    model<-list(y=y,Z=Z,H=H,T=T,R=R,Q=Q,a1=a1,P1=P1,P1inf=P1inf,u=e1$u,p=p,n=n,m=m,k=k,distribution=distribution,H_type=H_type,tolF=tolF,tol0=tol0)
    commontype<- intersect(class(e1),class(e2))
    if(length(commontype)==1){
        class(model)<-c("SSModel","customSSM")
    } else class(model)<-commontype
    
    invisible(model)
}