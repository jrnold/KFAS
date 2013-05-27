simDisturbance <- function(object, nsim = 1, antithetics = FALSE) {
    
    if (object$H_type == "Untransformed") 
        object <- transformSSM(object, type = "ldl")
    
    tv <- array(0, dim = 5)
    tv[1] <- dim(object$Z)[3] > 1
    tv[2] <- dim(object$H)[3] > 1
    tv[3] <- dim(object$T)[3] > 1
    tv[4] <- dim(object$R)[3] > 1
    tv[5] <- dim(object$Q)[3] > 1
    
    nsim2 <- 3 * nsim * antithetics + nsim
    ymiss <- array(is.na(object$y),dim=c(object$n,object$p))
    storage.mode(ymiss)<-"integer"
    
    if (object$H_type == "Augmented") {        
        
        etaplus <- array(0, c(object$r, object$n, nsim))
        aplus1 <- array(0, dim = c(object$m, nsim))
        
        x2 <- array(apply(object$Q, 3, diag) > object$tol0, c(object$r, (object$n - 1) * tv[5] + 1))
        x2 <- array(x2, c(object$r, object$n, nsim))
        dfeta <- sum(x2)/nsim
        
        nde <- which(diag(object$P1) > object$tol0)
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
            
            sims.out <- .Fortran("distsima", PACKAGE = "KFAS", NAOK = TRUE, ymiss = ymiss, tv=as.integer(tv), 
                    yt = array(object$y,dim=c(object$n,object$p)), Zt = object$Z, Tt = object$T, Rt = object$R, 
                    Qt = object$Q, a1 = object$a1, P1 = object$P1, object$P1inf, as.integer(nnd), as.integer(nsim),etaplus = etaplus, 
                    aplus1 = aplus1, p=object$p, n=object$n, m=object$m, r=object$r, info = as.integer(0), tolF = object$tolF, 
                    rankp = as.integer(sum(object$P1inf)), c = c2, etasim = array(0, c(object$r, object$n, nsim2)), object$tol0)
        } else {
            sims.out <- .Fortran("distsimnata", PACKAGE = "KFAS", NAOK = TRUE, ymiss = ymiss, tv=as.integer(tv), 
                    yt = array(object$y,dim=c(object$n,object$p)), Zt = object$Z, Tt = object$T, Rt = object$R, Qt = object$Q, 
                    a1 = object$a1, P1 = object$P1, object$P1inf, as.integer(nnd), as.integer(nsim), etaplus = etaplus, 
                    aplus1 = aplus1, p=object$p, n=object$n, m=object$m, r=object$r, info = as.integer(0), tolF = object$tolF, 
                    rankp = as.integer(sum(object$P1inf)), etasim = array(0, c(object$r, object$n, nsim2)), object$tol0)
        }
        sims <- list(eta = sims.out$etasim[, , 1:nsim, drop = FALSE])
    } else {
        
        
        
        epsplus <- array(0, c(object$p, object$n, nsim))
        etaplus <- array(0, c(object$r, object$n, nsim))
        aplus1 <- array(0, dim = c(object$m, nsim))
        
        x <- (array(apply(object$H, 3, diag) > object$tol0, c(object$p, object$n)))
        x <- array(x, c(object$p, object$n, nsim))
        dfeps <- sum(x)/nsim
        
        x2 <- array(apply(object$Q, 3, diag) > object$tol0, c(object$r, (object$n - 1) * tv[5] + 1))
        x2 <- array(x2, c(object$r, object$n, nsim))
        dfeta <- sum(x2)/nsim
        
        nde <- which(diag(object$P1) > object$tol0)
        nnd <- length(nde)
        dfu <- dfeps + dfeta + nnd
        
        u <- rnorm(dfu * nsim, mean = 0, sd = 1)
        
        if (dfeps > 0) 
            epsplus[x] <- u[1:(dfeps * nsim)]
        if (dfeta > 0) 
            etaplus[x2] <- u[(dfeps * nsim + 1):(dfeps * nsim + dfeta * nsim)]
        if (nnd > 0) 
            aplus1[nde, ] <- u[(dfeps * nsim + dfeta * nsim + 1):(dfu * nsim)]
        
        
        if (antithetics) {
            c2 <- numeric(nsim)
            for (i in 1:nsim) {
                u <- c(etaplus[, , i], epsplus[, , i], aplus1[, i])
                c2[i] <- t(u) %*% c(u)
            }
            q <- pchisq(c2, df = dfu)
            c2 <- sqrt(qchisq(1 - q, dfu)/c2)
            sims.out <- .Fortran("distsimd", PACKAGE = "KFAS", NAOK = TRUE, ymiss = ymiss, tv=as.integer(tv), 
                    yt = array(object$y,dim=c(object$n,object$p)), Zt = object$Z, H = object$H, Tt = object$T, 
                    Rt = object$R, Qt = object$Q, a1 = object$a1, P1 = object$P1, object$P1inf, as.integer(nnd), 
                    as.integer(nsim), epsplus = epsplus, etaplus = etaplus, aplus1 = aplus1, object$p, object$n, object$m, 
                    object$r, info = as.integer(0), tolF = object$tolF, rankp = as.integer(sum(object$P1inf)), 
                    c = c2, epssim = array(0, c(object$p, object$n, nsim2)), etasim = array(0, c(object$r, object$n, nsim2)), object$tol0)
        } else {
            sims.out <- .Fortran("distsimnatd", PACKAGE = "KFAS", NAOK = TRUE, ymiss = ymiss, tv=as.integer(tv), 
                    yt = array(object$y,dim=c(object$n,object$p)), Zt = object$Z, H = object$H, Tt = object$T, Rt = object$R, 
                    Qt = object$Q, a1 = object$a1, P1 = object$P1, object$P1inf, as.integer(nnd), as.integer(nsim), 
                    epsplus = epsplus, etaplus = etaplus, aplus1 = aplus1, p=object$p, n=object$n, m=object$m, r=object$r, 
                    info = as.integer(0), tolF = object$tolF, rankp = as.integer(sum(object$P1inf)), 
                    epssim = array(0, c(object$p, object$n, nsim2)), etasim = array(0, c(object$r, object$n, nsim2)), object$tol0)
            
        }
        sims <- list(eps = sims.out$epssim, eta = sims.out$etasim)
    }
    
    if (sims.out$info != 0) {
        if (sims.out$info == 1) 
            stop("Couldn't compute LDL decomposition of Ht!")
        if (sims.out$info == 2) 
            stop("Couldn't compute LDL decomposition of Qt!")
        if (sims.out$info == 3) 
            stop("Couldn't compute LDL decomposition of P1!")
    }
    
    sims
} 
