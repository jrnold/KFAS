library(KFAS)
    
kk <- c(396, 399, 403, 434, 340, 519, 558, 566, 591, 566, 574, 646, 644,
        646, 665, 693, 773, 834, 910,1035, 1002, 1161, 1056, 1097, 1094, 1042,
        1194, 1316, 1246, 1503, 1428, 1477, 1490, 1465, 1560, 1860, 2008, 2020,
        2167) 

vlk <- c(4623785, 4606307, 4612124, 4639657, 4666081, 4690574, 4711440,
        4725664, 4738902, 4752528, 4764690, 4779535, 4799964, 4826933, 4855787,
        4881803, 4902206, 4918154, 4932123, 4946481, 4964371, 4986431, 5013740,
        5041992, 5066447, 5088333, 5107790, 5124573, 5139835, 5153498, 5165474,
        5176209, 5188008, 5200598, 5213014, 5228172, 5246096, 5266268, 5288720)


n<-39
Zt<-array(c(1,0),c(1,2))
Tt<-diag(c(1,1))
Tt[1,2]<-1
Rt <- diag(2)

model<-SSModel(y=kk,u=vlk,Z=Zt,T=Tt,Q=diag(NA,2),R=Rt,distribution="Poisson")
likfn<-function(pars)
{
    diag(model$Q[,,1])<-exp(pars)
    logLik(model)
}

opt <- optim(par=c(1, 1), fn=likfn,  method="BFGS",  control=list(fnscale=-1, trace=1, REPORT=1))

diag(model$Q[,,1])<-exp(opt$par)

out<-KFS(model,nsim=1000)
mod<-model
mod$y[]<-NA
mod$a1[]<-out$alphahat[,1]
mod$P1inf[]<-0
mod$P1[]<-0
isample<-importanceSSM(mod,nsim=1000)

yy<-t(t(vlk*exp(isample$s[1,,]))*isample$w)
ts.plot(cbind(yy,model$y),col=c(rep(2,4000),1))
###


# Example of model building function

# Number of alcohol related deaths in Finland in a certain age cohort.
# Local level model with a drift

y<-c(0,6,2,1,2,1,2,3,1,0,4,2,1,2,2,8,4,2,1,1,2,0,3,3,1,1,3,3,4,3,3,0,6,5,3,5,2,3,3)

u<-c(0.042580,0.043765,0.045320,0.047410,0.049765,0.051620,0.053130,0.055340,
0.058575,0.061610,0.065000,0.069005,0.072880,0.077430,0.082835,0.089210,
0.094010,0.098685,0.104930,0.110690,0.116645,0.121965,0.127640,0.134900,
0.142085,0.150540,0.159485,0.167140,0.174490,0.179785,0.183535,0.186870,
0.189335,0.192170,0.196285,0.199770,0.208725,0.223460,0.237530)


model<-structSSM(y=y,trend="slope",Q.slope=0,u=u,distribution="Poisson")


fit<-fitSSM(model=model,inits=-2)
model<-fit$model

amod<-approxSSM(model)

out<-KFS(model,smoothing="state",nsim=1000)

ts.plot(amod$theta)
lines(signal(out)$signal,col=3)
# approximating model seems to overestimate the signal

plot(amod) # u*exp(amod$theta), approximation
lines(out$model$u*exp(signal(out)$signal),col=3) # this is biased way
lines(out$yhat,col=4) # This is correct way
# approximating model gives very close values

# It looks like the model gives too smooth signal, let's add additional noise term to scope with overdispersion:

fn<-function(par, y, u){
    
    Z <- matrix(c(1,0,1), 1,3)
    T <- matrix(c(1,0,0, 1,1,0, 0,0,0), 3,3)
    Q <- diag(c(exp(par[1]),0,exp(par[2])))
    P1inf<-diag(c(1,1,0))
    P1<-diag(c(0,0,exp(par[2])))
    SSModel(y = y, Z = Z, T = T, Q = Q, u = u, P1=P1,P1inf=P1inf,dist = "Poisson")
    
}

fit<-fitSSM(modFun=fn,inits=c(-2,-2),y=y,u=u)
fit$model$Q
out<-KFS(fit$model,smoothing="state",nsim=1000)
lines(out$yhat,col=5) 

###


data(boat)

model<-structSSM(y=boat,distribution="Binomial")

fit<-fitSSM(model,inits=log(0.5))

exp(fit$opt$p)

model<-fit$model

model$Q[]<-0.521
out<-KFS(model)
ts.plot(out$yhat,col=1,ylim=c(0,1))
points(model$y,pch=19)