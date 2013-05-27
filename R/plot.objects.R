#' Plot Ouput of Kalman Filter and Smoother
#' @S3method plot KFS
#' @method plot KFS
#' @param x object of class \code{KFS}
#' @param y Ignored.
#' @param type Draw \code{"signals"} together with observations, or figures of \code{"state"} components. 
#' Default is \code{"signal"}. If smoothed values of states are available, those are used. Otherwise filtered values are used.
#' @param ... Graphical parameters for function plot.ts.

plot.KFS<-function(x, y=NULL, type="signal",...){
    type<-match.arg(type,choices=c("signal","state"))
    if(type=="signal"){       
        Series<-ts.union(x$model$y,signal(x)$signal)        
        plot.type<-"single"        
    } else {
        if(!is.null(x$alphahat)){                            
            Series<-ts(t(x$alphahat),start=start(x$model$y))            
        } else {            
            Series<-ts(t(x$a)[-(x$model$n+1),],start=start(x$model$y))
        }
        plot.type<-"multiple"
    }   
    
    plot.ts(Series,plot.type=plot.type,...) 
    
    
}
#' Plot Approximating Gaussian State Space Model
#' @S3method plot approxSSM
#' @method plot approxSSM
#' @param x approxSSM object
#' @param y Ignored.
#' @param legend.position Position of the legend as a keyword. See function legend for details.
#' @param ... Graphical parameters for function plot.ts.

plot.approxSSM<-function(x, y=NULL,legend.position="bottomright",...){    
    
    if(x$original.distribution=="Poisson"){        
        mean<-x$u*exp(x$thetahat)
        meantext<-"Mean u*exp(thetahat)"
    } else {
        mean<-x$u*exp(x$thetahat)/(1+exp(x$thetahat))
        meantext<-"Mean u*exp(thetahat)/(1+exp(thetahat))"
    }
    Series<-ts.union(x$original.y,mean) 
    plot.ts(Series,plot.type="single", col=1:2,...) 
    legend(x=legend.position,y=NULL,legend=c("Observations y",meantext), col=1:2, lty=1)    
    
}
