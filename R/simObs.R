simObs <- function(object, nsim = 1, antithetics = FALSE) {
    
    
    if (sum(is.na(object$y)) == 0) 
        stop("There is no missing observations, nothing to simulate.")
    
    tv <- array(0, dim = 5)
    tv[1] <- dim(object$Z)[3] > 1
    tv[2] <- dim(object$H)[3] > 1
    tv[3] <- dim(object$T)[3] > 1
    tv[4] <- dim(object$R)[3] > 1
    tv[5] <- dim(object$Q)[3] > 1
    ymiss <- array(is.na(object$y),dim=c(object$n,object$p))
    storage.mode(ymiss)<-"integer"    
    
    if (identical(object$H_type, "LDL decomposed"))
        warning("Because model is already LDL-transformed, simulated series are are also in transformed scale.")

    if (! identical(object$H_type, "Augmented"))
        object <- transformSSM(object, type = "augment")
    
    
    nsim2 <- 3 * nsim * antithetics + nsim
    etaplus <- array(0, c(object$k, object$n, nsim))
    aplus1 <- array(0, dim = c(object$m, nsim))
    
    x2 <- array(apply(object$Q, 3, diag) > object$tol0, c(object$k, (object$n - 1) * tv[5] + 1))
    x2 <- array(x2, c(object$k, object$n, nsim))
    dfeta <- sum(x2)/nsim
    
    nde <- which(diag(object$P1) > object$tol0)
    nd <- which(diag(object$P1inf) == 0)
    nnd <- length(nde)
    dfu <- dfeta + nnd
    
    u <- rnorm(dfu * nsim, mean = 0, sd = 1)
    
    if (dfeta > 0) 
        etaplus[x2] <- u[1:(dfeta * nsim)]
    if (nnd > 0) 
        aplus1[nde, ] <- u[(dfeta * nsim + 1):(dfu * nsim)]
    
    if (antithetics) {
        c2 <- numeric(nsim)
        for (i in 1:nsim) {
            u <- c(etaplus[, , i], aplus1[, i])
            c2[i] <- t(u) %*% c(u)
        }
        q <- pchisq(c2, df = dfu)
        c2 <- sqrt(qchisq(1 - q, dfu)/c2)
        
        sims.out <- .Fortran("obssim", PACKAGE = "KFAS", NAOK = TRUE, ymiss = ymiss, 
            tv=as.integer(tv), yt = array(object$y,dim=c(object$n,object$p)), Zt = object$Z, 
            Tt = object$T, Rt = object$R, Qt = object$Q, a1 = object$a1, P1 = object$P1, 
            object$P1inf, as.integer(nnd), as.integer(nsim), etaplus = etaplus, aplus1 = aplus1, 
            p=object$p, n=object$n, m=object$m, r=object$k, info = integer(1), 
            tolF = object$tolF, rankp = as.integer(sum(object$P1inf)), c = c2, 
            ysim = array(object$y, c(object$n, object$p, nsim2)), object$tol0,as.integer(nd),as.integer(length(nd)))
    } else {
        sims.out <- .Fortran("obssimnat", PACKAGE = "KFAS", NAOK = TRUE, ymiss = ymiss, 
            tv=as.integer(tv), yt = array(object$y,dim=c(object$n,object$p)), Zt = object$Z, 
            Tt = object$T, Rt = object$R, Qt = object$Q, a1 = object$a1, P1 = object$P1, 
            object$P1inf, as.integer(nnd), as.integer(nsim), etaplus = etaplus, aplus1 = aplus1, 
            p=object$p, n=object$n, m=object$m, r=object$k, info = integer(1), 
            tolF = object$tolF, rankp = as.integer(sum(object$P1inf)), 
            ysim = array(object$y, c(object$n, object$p, nsim2)), object$tol0,as.integer(nd),as.integer(length(nd)))
        
    }
    if (sims.out$info != 0) {
        if (sims.out$info == 2) 
            stop("Couldn't compute LDL decomposition of Qt! Try changing the tol0 parameter of the model.")
        if (sims.out$info == 3) 
            stop("Couldn't compute LDL decomposition of P1! Try changing the tol0 parameter of the model.")
    }
    sims.out$ysim
} 
