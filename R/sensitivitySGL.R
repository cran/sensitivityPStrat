.calc.time.seq.list <- function(time.points, times) {
  calc.time.seq <- function(time.point, times) c(times <= time.point, FALSE)
  
  q.seq.list <- lapply(time.points, times, FUN=calc.time.seq)
  names(q.seq.list) <- time.points

  return(q.seq.list)
}

.calc.time.seq.index <- function(time.points, times) {
  calc.time.seq <- function(time.point, times) sum(times <= time.point)
  
  q.index <- unlist(lapply(time.points, times, FUN=calc.time.seq))
  names(q.index) <- time.points

  return(q.index)
}

.calcSGLF1 <- function(time.points, KM1) {
  q.index <- .calc.time.seq.index(time.points, KM1$t)

  ans <- KM1[ifelse(q.index == 0, NA, q.index),]
  ans[q.index == 0,] <- 0
  
  return(ans)
}

.calc.beta.tplus.list <- function(beta, t0, tau) {
  tplus0 <- c(pmin(t0, tau), tau)

  calc.beta.tplus <- function(beta,time)
    if(is.finite(beta)) return(beta*time) else return(beta)
  return(lapply(beta, time=tplus0, FUN=calc.beta.tplus))
}

.calcSGLCoeffCommon <- function(beta, KM0, tau, time.points) {
  q.list <- .calc.time.seq.list(time.points, KM0$t)
  q.index <- unlist(lapply(q.list, FUN=sum))

  beta.tplus <- .calc.beta.tplus.list(beta, KM0$t, tau)

  KMAns <- KM0[ifelse(q.index == 0L, NA, q.index),]
  KMAns[q.index == 0L,] <- 0

  return(list(q.list=q.list, q.index=q.index, beta.tplus=beta.tplus,
              KMAns=KMAns, i=seq_along(beta.tplus)))
}

.calcSGL.coeff.smp <- function(beta, KM0, dF0, tau, time.points, RR) {
  com <- .calcSGLCoeffCommon(beta, KM0, tau, time.points)
  
  return(mapply(FUN=.calcSGLBetaCoeffBasic,
                i=com$i, beta=beta, beta.tplus=com$beta.tplus,
                MoreArgs=list(q.index=com$q.index, KMAns=com$KMAns, dF0=dF0,
                  RR=RR),
                SIMPLIFY=FALSE, USE.NAMES=FALSE))
}

.calcSGL.coeff.adv <- function(beta, KM0, dF0, p0, t0, n0, N0, n1, N1, tau,
                           time.points, RR) {

  com <- .calcSGLCoeffCommon(beta, KM0, tau, time.points)

  coeffs <- mapply(FUN=.calcSGLBetaCoeffAdv,
                   i=com$i, beta=beta, beta.tplus=com$beta.tplus,
                   MoreArgs=list(q.list=com$q.list, q.index=com$q.index,
                     KMAns=com$KMAns, dF0=dF0, p0=p0, n0=n0, N0=N0,
                     n1=n1, N1=N1, RR=RR, len.F=length(KM0$t)),
                   USE.NAMES=FALSE, SIMPLIFY=FALSE)

  return(coeffs)
}

.calcSGL.g <- function(q.index, q.seq, w, w.dF, w.1mw.dF,
                   sum.w.dF, sum.w.1mw.dF, diff.w, Fas) {

  sum.w.dF.qseq <- sum(w.dF[q.seq])
  sum.w.dF.not.qseq <- sum(w.dF[!q.seq])

  dg <- c(0,
          ## dgdalpha
          (sum.w.dF * sum(w.1mw.dF[q.seq]) - sum.w.1mw.dF*sum.w.dF.qseq),
          ## dgdF where j < l
          diff.w[replace(q.seq, q.index, FALSE)] * sum.w.dF.not.qseq,
          ## dgdF where j == l
          (sum.w.dF.not.qseq*w[q.index] + sum.w.dF.qseq*w[q.index+1]),
          ## dgdF where j > l
          sum.w.dF.qseq*diff.w[!q.seq[-length(q.seq)]]) / (sum.w.dF * sum.w.dF)

  return(dg)
}

.calcSGLPosInfFas <- function(KMAns, RR)
  pmin(KMAns$Fas/RR, 1L)

.calcSGLNegInfFas <- function(KMAns, RR)
  pmax((KMAns$Fas - 1L)/RR + 1L, 0)

.calcSGLPosInfFasVar <- function(KMAns, RR, n0, N0, n1, N1)
  (KMAns$Fas.var/RR)^2 + (KMAns$Fas/RR)^2*(1L/n0 - 1L/N0 + 1L/n1 - 1L/N1)

.calcSGLNegInfFasVar <- function(KMAns, RR, n0, N0, n1, N1)
  (KMAns$Fas.var/RR)^2 + ((1 - KMAns$Fas)/RR)^2*(1/n0 - 1/N0 + 1/n1 - 1/N1)

.calcSGLInfCoeffAdv <- function(beta, KMAns, RR, n0, N0, n1, N1) {

  if(beta < 0L) {
    Fas <- .calcSGLNegInfFas(KMAns, RR)
    Fas.var <- .calcSGLNegInfFasVar(KMAns, RR, n0, N0, n1, N1)
  } else {
    Fas <- .calcSGLPosInfFas(KMAns, RR)
    Fas.var <- .calcSGLPosInfFasVar(KMAns, RR, n0, N0, n1, N1)
  }

  return(list(Fas=Fas, Fas.var=Fas.var))
}

.calcSGLInfCoeffBasic <- function(beta, KMAns, RR) {
  if(beta < 0L) {
    Fas <- .calcSGLNegInfFas(KMAns, RR)
  } else {
    Fas <- .calcSGLPosInfFas(KMAns, RR)
  }

  return(list(Fas=Fas))
}

.calc.w <- function(alpha, beta.tplus) 
  1/(1L+exp(-alpha - beta.tplus))

.calc.alphahat <- function(beta.tplus, dF, B) {
  ee <- function(a, beta.tplus, dF, B) {
    (sum(.calc.w(alpha=a, beta.tplus=beta.tplus)*dF) - B)^2
  }

  alphahat <- optimize(ee, c(-100,100), beta.tplus=beta.tplus, dF=dF,
                       B=B)$minimum

  if(alphahat > 90 || alphahat < -90) {
    warning("optimize overflow alphahat value invalid")
  }

  alphahat
}

.calcSGLBetaCoeffBasic <- function(i, beta, beta.tplus, q.index, KMAns, dF0,
                                  RR) {
  if(is.infinite(beta)) {
    return(c(list(i=i, alphahat=NA),
             .calcSGLInfCoeffBasic(beta, KMAns, RR)))
  }

  alphahat <- .calc.alphahat(beta.tplus, dF0, RR)

  if(all.equal(beta, 0) == TRUE)
    return(list(i=i, alphahat=alphahat, Fas=KMAns$Fas))

  w <- .calc.w(alpha=alphahat, beta.tplus=beta.tplus)

  w.dF0 <- w*dF0
  F0ai <- cumsum(w.dF0)/RR

  Fas <- ifelse(q.index == 0L, 0L, F0ai[q.index])

  return(list(i=i, alphahat=alphahat, Fas=Fas))
}

.calcSGLBetaCoeffAdv <- function(i, beta, beta.tplus, q.list, q.index, KMAns, dF0,
                             p0, n0, N0, n1, N1, RR, len.F) {
  if(is.infinite(beta)) {
    return(c(list(i=i, alphahat=NA),
             .calcSGLInfCoeffAdv(beta, KMAns, RR, n0, N0, n1, N1)))
  }
  
  alphahat <- .calc.alphahat(beta.tplus, dF0, RR)

  if(all.equal(beta, 0) == TRUE) {
    Fas <- KMAns$Fas

    w <- .calc.w(alpha=alphahat, beta.tplus=0)

    A1 <- -N1*p0*w*(1-w)
    B1 <- 0

    dg.nrow <- len.F + 2L
    dg <- matrix(0, ncol=length(q.index), nrow=dg.nrow,
                 dimnames=list(NULL, names(q.list)))
    
    dg[seq.int(from=0L, along=q.index)*dg.nrow + q.index + 2] <- 1
    dg <- as.data.frame(dg)
    
    return(list(i=i, alphahat=alphahat, A1=A1, B1=B1, Fas=Fas, Fas.var=NULL,
                dg=dg))
  }
  
    
  w <- .calc.w(alpha=alphahat, beta.tplus=beta.tplus)

  w.dF0 <- w*dF0
  F0ai <- cumsum(w.dF0)/RR
  Fas <- ifelse(q.index == 0, 0, F0ai[q.index])
  
  w.1mw.dF0 <- w.dF0*(1-w)
  sum.w.dF0 <- sum(w.dF0)
  sum.w.1mw.dF0 <- sum(w.1mw.dF0)
  diff.w <- diff(w)
  
  A1 <- sum(-N1 * w.1mw.dF0 * p0)
  B1 <- -N1*p0*-diff.w

  dg <- mapply(FUN=.calcSGL.g, q.index=q.index, q.seq=q.list,
                 MoreArgs=list(w, w.dF0, w.1mw.dF0,
                 sum.w.dF0, sum.w.1mw.dF0, diff.w, F0ai),
                 SIMPLIFY=FALSE)

  return(list(i=i, A1=A1, B1=B1, alphahat=alphahat, Fas=Fas, Fas.var=NULL,
              dg=dg))
}

.makeBootstrapLenIndx <- function(s, indx.seq, N)
  sample(indx.seq, N, replace=TRUE)

sensitivitySGL <- function(z, s, d, y, beta, tau, time.points,
                           selection, trigger, groupings,
                           empty.principal.stratum,
                           ci=0.95, ci.method=c("analytic", "bootstrap"),
                           na.rm=FALSE, N.boot=100L, oneSidedTest=FALSE,
                           twoSidedTest=TRUE, verbose=getOption("verbose")) {
  if(!require(survival))
    stop("require's the survival package to function")

  ## z - group that subject belongs to
  ## s - subject met selection cirteria
  ## d - subject had event
  ## y - time until event ocurred

  if(!missing(ci.method) && is.null(ci.method))
    isSlaveMode <- TRUE
  else
    isSlaveMode <- FALSE

  if(!isSlaveMode) {
    ## Not running a boot strap mode
    ## Run error checks on variables.
    ErrMsg <- c(.CheckEmptyPrincipalStratum(empty.principal.stratum),
                .CheckSelection(selection, s, empty.principal.stratum),
                .CheckGroupings(groupings),
                .CheckTrigger(trigger, d),
                .CheckLength(z=z, s=s, y=y),
                .CheckZ(z, groupings, na.rm),
                .CheckS(s, empty.principal.stratum, na.rm),
                .CheckY(y, s, selection),
                .CheckD(d=d, s=s, selection=selection))

    if(length(ErrMsg) > 0L)
      stop(paste(ErrMsg, collapse="\n  "))


    if(na.rm == TRUE) {
      naIndex <- !(is.na(s) | is.na(z) | (s & (is.na(d) | is.na(y))))

      z <- z[naIndex]
      s <- s[naIndex]
      d <- d[naIndex]
      y <- y[naIndex]
    }

    GroupReverse <- FALSE
    if(empty.principal.stratum[1L] == selection) {
      z <- ifelse(z == groupings[1L], TRUE, ifelse(z == groupings[2L], FALSE, NA))
      GroupReverse <- TRUE
    } else if(empty.principal.stratum[2L] == selection)
      z <- ifelse(z == groupings[2L], TRUE, ifelse(z == groupings[1L], FALSE, NA))
    s <- s == selection
    d <- d == trigger
  } else {
    GroupReverse <- ci
  }
  
  ci.method <- sort(unique(match.arg(ci.method, several.ok=TRUE)))
  
  ## N  - number subjects
  ## N0 - number of subjects in group 0
  ## N1 - number of subjects in group 1
  N <-length(z)
  N1 <- sum(z)
  N0 <- N-N1

  z0.s1 <- !z & s
  z1.s1 <- z & s
  
  ## n0 - number of subjects in group 0 that were selected 
  ## n1 - number of subjects in group 1 that were selected
  n0 <- sum(z0.s1)
  n1 <- sum(z1.s1)

  ## p0 - probiblity that subject in group 0 was selected
  ## p1 - probablity that subject in group 1 was selected
  p0 <- n0/N0
  p1 <- n1/N1

  RR <- p1/p0
  VE <- 1L - RR

  ## summary survfit of length of time til event for group 0 and group 1.
  KM0 <- with(summary(survfit(Surv(y[z0.s1],d[z0.s1])~1L)),
              data.frame(t=time, Fas=1L - surv, Fas.var=std.err*std.err))

  KM1 <- with(summary(survfit(Surv(y[z1.s1],d[z1.s1])~1L)), 
              data.frame(t=time, Fas=1L - surv, Fas.var=std.err*std.err))

  len.t0 <- length(KM0$t)
  
  ## total length of phi vector
  len.total <- 2L + len.t0
  
  dF0 <- diff(c(0L,KM0$Fas,1L))

  doAnalyticCi <- ('analytic' %in% ci.method && !is.null(ci))
  doBootStrapCi <- ('bootstrap' %in% ci.method && !is.null(ci))

  betaOrig <- beta
  beta <- unique(sort(beta))

  timePointsOrig <- time.points
  time.points <- unique(sort(time.points))

  if(doAnalyticCi) {
    coeffs0 <- .calcSGL.coeff.adv(beta=beta, KM0=KM0, dF0=dF0, p0=p0,
                              n0=n0, N0=N0, n1=n1, N1=N1, tau=tau,
                                 time.points=time.points, RR=RR)
  } else {
    coeffs0 <- .calcSGL.coeff.smp(beta=beta, KM0=KM0, dF0=dF0, tau=tau,
                                   time.points=time.points, RR=RR)
  }

  coeff1 <- .calcSGLF1(time.points=time.points, KM1)
  
  SCE.dim <- c(length(time.points), length(beta))
  SCE.length <- prod(SCE.dim)
  SCE.dimnames <- list(time.points=as.character(time.points),
                       beta=as.character(beta))
  SCE <- array(numeric(SCE.length), dim=SCE.dim, dimnames=SCE.dimnames)
  
  for(coeff0 in coeffs0) {
    SCE[,coeff0$i] <- coeff0$Fas - coeff1$Fas
  }
  
  if(isSlaveMode) return(list(SCE=SCE))
  
  SCE.var.dim <- c(SCE.dim, length(ci.method))
  SCE.var.dimnames <- c(SCE.dimnames, list(ci.method=ci.method))
  
  SCE.var <- array(numeric(0),
                    dim=SCE.var.dim,
                    dimnames=SCE.var.dimnames)

  if(doAnalyticCi) {
    Gamma <- Omega <- matrix(0, ncol=len.total, nrow=len.total)

    b.col1 <- seq(from=3, length=len.t0)

    ##  N
    ## ====
    ## \
    ##  >    (1 - Z ) (S  - p0)
    ## /           i    i
    ## ====
    ## i = 1
    ##
    ## this can be broken down into the sum of all values where (1 - Z[i]) is 1
    ## (1 to N0) plus the sum of all remaining values where (1 - Z[i]) is 0.
    ## The value of the sum
    ## of all values where (1 - Z[i]) is 0 will alway be zero.
    ##
    ##  N                              N0
    ## ====                           ====
    ## \                              \
    ##  >         (1 - 1) (S  - p0) +  >    (1 - 0) (S  - p0)
    ## /                    i         /              i
    ## ====                           ====
    ## i = N0 + 1                     i = 1
    ##
    ## Simplifying to
    ##
    ##  N0
    ## ====
    ## \
    ##  >    (S  - p0)
    ## /       i
    ## ====
    ## i = 1
    ##
    ## This can then be further broken down in the sum of all remaining values
    ## (from 1 to N0) where S[i] is 1 (from 1 to n0); plus the sum of all
    ## remaining values where S[i] is 0 (from n0+1 to N0).
    ##
    ##  n0               N0
    ## ====             ====
    ## \                \          
    ##  >    (1 - p0) +  >    (0 - p0)
    ## /                /          
    ## ====             ====
    ## i = 1            i = n0 + 1
    ##
    ## Simplifying to
    ##
    ## n0 * (1 - p0) - p0 * (N0 - n0)
    ##
    ## p0*N0 - n0
    ##
    Omega[1,1] <- p0*(N0 - n0)

    Omega[2,2] <- n1 - n1*p1
    ## Calcuate the V for z1.s1
    V <- calc.v(d[z0.s1], y[z0.s1])
    VmF0 <- V - KM0$Fas

    Omega[1, b.col1] <- rowSums(VmF0 * (1-p0))
    .sumCrossUpperTri(Omega[b.col1,b.col1]) <- VmF0
    
    ## Fold Omega
    Omega <- .foldUpperTri(Omega)

    ## populate Gamma matrix
    diag(Gamma) <- rep(c(-N0, -N1, -n0), times=c(1,1,len.t0))

    Gamma[1,2] <- -N1*RR

    for(coeff0 in coeffs0) {
      if(is.infinite(beta[coeff0$i])) {
        SCE.var[,coeff0$i,'analytic'] <- coeff0$Fas.var
      } else {
        Gamma[2,2] <- coeff0$A1
        Gamma[b.col1,2] <- coeff0$B1

        IGamma <- solve(Gamma/N)
        vartheta <- crossprod(IGamma, Omega) %*% IGamma / N

        SCE.var[,coeff0$i,'analytic'] <- sapply(coeff0$dg, function(dg) {
          (dg %*% vartheta %*% dg) }) / N + coeff1$Fas.var
      }
    }
  }
  
  if(twoSidedTest) {
    ci.probs <- c(ifelse(ci < 0.5, ci, 0), ifelse(ci < 0.5, 1, ci)) +
      rep(ifelse(ci < 0.5, -ci/2, (1-ci)/2), times=length(ci))
  } else {
    ci.probs <- NULL
  }

  if(oneSidedTest) {
    ci.probs <- c(ci.probs, ci)
  }

  ci.probs <- unique(ci.probs)
  ci.probsLen <- length(ci.probs)

  SCE.ci.dim <- c(ci.probsLen, SCE.var.dim)
  SCE.ci.dimnames <- c(list(ci.probs=as.character(ci.probs)), SCE.var.dimnames)
  
  SCE.ci <- array(numeric(0),
                  dim=SCE.ci.dim,
                  dimnames=SCE.ci.dimnames)
  
  if(doAnalyticCi) {
    SCE.rep <- rep.int(ci.probsLen, times=SCE.length)

    SCE.ci[,,,'analytic'] <- array(rep.int(SCE, times=SCE.rep) + 
                                   rep.int(qnorm(ci.probs), times=SCE.length) * 
                                   rep.int(sqrt(SCE.var[,,'analytic']),
                                           times=SCE.rep),
                                   dim=c(ci.probsLen, SCE.dim))
  }
  
  if(doBootStrapCi) {
    
    z.seq <- seq_len(N)
    current.fun <- sys.function()

    bootCalc <- function(i, z.seq, nVal, beta, tau, time.points,
                         current.fun, verbose) {
      samp <- .makeBootstrapLenIndx(s, indx.seq=z.seq, N=nVal)
      ans <- as.vector(current.fun(z=z[samp], s=s[samp], d=d[samp],
                                   y=y[samp],
                                   beta=beta,
                                   tau=tau,
                                   selection=TRUE,
                                   groupings=c(FALSE,TRUE),
                                   empty.principal.stratum=c(FALSE,TRUE),
                                   trigger=TRUE,
                                   time.points=time.points,
                                   ci.method=NULL,
                                   ci=GroupReverse)$SCE)
      if(verbose) cat(".")
      return(ans)
    }

    vals <- apply(do.call(rbind, lapply(integer(N.boot), FUN=bootCalc,
                                        z.seq=z.seq, nVal=N,
                                        beta=beta,
                                        tau=tau,
                                        time.points=time.points,
                                        current.fun=current.fun,
                                        verbose=verbose)),
                  2,
                  FUN=function(x) return(c(var(x), quantile(x, probs=ci.probs))))
    if(verbose) cat("\n")

    SCE.var.boot = vals[1,,drop=FALSE]
    SCE.ci.boot = vals[-1,,drop=FALSE]
    
    dim(SCE.var.boot) <- SCE.dim
    dim(SCE.ci.boot) <- c(ci.probsLen, SCE.dim)

    SCE.var[,,'bootstrap'] <- SCE.var.boot

    str(SCE.ci.boot)
    str(SCE.ci)
    SCE.ci[,,,'bootstrap'] <- SCE.ci.boot
  }

  tpIndex <- match(time.points, timePointsOrig)
  betaIndex <- match(beta, betaOrig)
  ans <- list(SCE=SCE[tpIndex,betaIndex,drop=FALSE],
              SCE.var=SCE.var[tpIndex, betaIndex, ,drop=FALSE],
              SCE.ci=SCE.ci[,tpIndex, betaIndex, , drop=FALSE],
              beta=betaOrig)
  
  class(ans) <- "plot.sensitivity2.5d"

  return(ans)
}
