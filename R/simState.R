simState <- function(object, nsim = 1, antithetics = FALSE) {
    
    if (object$H_type == "Untransformed") 
        object <- transformSSM(object = object, type = "ldl")
    
    tv <- array(0, dim = 5)
    tv[1] <- dim(object$Z)[3] > 1
    tv[2] <- dim(object$H)[3] > 1
    tv[3] <- dim(object$T)[3] > 1
    tv[4] <- dim(object$R)[3] > 1
    tv[5] <- dim(object$Q)[3] > 1

    ymiss <- array(is.na(object$y),dim=c(object$n,object$p))
    storage.mode(ymiss)<-"integer"
    
    nsim2 <- 3 * nsim * antithetics + nsim
    
    
    epsplus <- array(0, c(object$p, object$n, nsim))
    etaplus <- array(0, c(object$k, object$n, nsim))
    aplus1 <- array(0, dim = c(object$m, nsim))
     
    
    x <- (array(apply(object$H, 3, diag) > object$tol0, c(object$p, object$n)) & (!t(ymiss)))
    x <- array(x, c(object$p, object$n, nsim))
    dfeps <- sum(x)/nsim
    
    x2 <- array(apply(object$Q, 3, diag) > object$tol0, c(object$k, (object$n - 1) * tv[5] + 1))
    x2 <- array(x2, c(object$k, object$n, nsim))
    dfeta <- sum(x2)/nsim
    
    nde <- which(diag(object$P1) > object$tol0)
    nd <- which(diag(object$P1inf) == 0)
    nnd <- length(nde)
    dfu <- dfeps + dfeta + nnd
    u <- rnorm(dfu * nsim, mean = 0, sd = 1)
    
    if (dfeps > 0) 
        epsplus[x] <- u[1:(dfeps * nsim)]
    if (dfeta > 0) 
        etaplus[x2] <- u[(dfeps * nsim + 1):(dfeps * nsim + dfeta * nsim)]
    if (nnd > 0) 
        aplus1[nde, ] <- u[(dfeps * nsim + dfeta * nsim + 1):(dfu * nsim)]
    
    c2 <- numeric(nsim)
    
    if (antithetics) {
        for (i in 1:nsim) {
            u <- c(etaplus[, , i], epsplus[, , i], aplus1[, i])
            c2[i] <- t(u) %*% c(u)
        }
        q <- pchisq(c2, df = dfu)
        c2 <- sqrt(qchisq(1 - q, dfu)/c2)
        
        sims.out <- .Fortran("alphasim", PACKAGE = "KFAS", NAOK = TRUE, ymiss = ymiss, 
            tv=as.integer(tv), yt = array(object$y,dim=c(object$n,object$p)), Zt = object$Z, H = object$H, Tt = object$T, Rt = object$R, 
            Qt = object$Q, a1 = object$a1, P1 = object$P1, object$P1inf, as.integer(nnd), 
            as.integer(nsim), epsplus = epsplus, etaplus = etaplus, aplus1 = aplus1, 
            object$p, object$n, object$m, object$k, info = as.integer(0), tolF = object$tolF, rankp = as.integer(sum(object$P1inf)), 
            c = c2, alphasim = array(0, c(object$m, object$n, nsim2)), object$tol0,as.integer(nd),as.integer(length(nd)))
    } else {
        sims.out <- .Fortran("alphasimnat", PACKAGE = "KFAS", NAOK = TRUE, ymiss = ymiss, 
            tv=as.integer(tv), yt = array(object$y,dim=c(object$n,object$p)), Zt = object$Z, H = object$H, Tt = object$T, Rt = object$R, 
            Qt = object$Q, a1 = object$a1, P1 = object$P1, object$P1inf, as.integer(nnd), 
            as.integer(nsim), epsplus = epsplus, etaplus = etaplus, aplus1 = aplus1, 
            object$p, object$n, object$m, object$k, info = as.integer(0), tolF = object$tolF, rankp = as.integer(sum(object$P1inf)), 
            alphasim = array(0, c(object$m, object$n, nsim2)), object$tol0,as.integer(nd),as.integer(length(nd)))
    }
    if (sims.out$info != 0) {
        if (sims.out$info == 1) 
            stop("Couldn't compute LDL decomposition of Ht!")
        if (sims.out$info == 2) 
            stop("Couldn't compute LDL decomposition of Qt!")
        if (sims.out$info == 3) 
            stop("Couldn't compute LDL decomposition of P1!")
    }
    
    sims.out$alphasim
} 
