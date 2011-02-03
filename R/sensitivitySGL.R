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
  return(list(beta.tplus=lapply(beta, time=tplus0, FUN=calc.beta.tplus),
              tplus=tplus0))
}

.calcSGLCoeffCommon <- function(beta, KM0, tau, time.points) {
  q.list <- .calc.time.seq.list(time.points, KM0$t)
  q.index <- unlist(lapply(q.list, FUN=sum))

  tvals <- .calc.beta.tplus.list(beta, KM0$t, tau)

  KMAns <- KM0[ifelse(q.index == 0L, NA, q.index),]
  KMAns[q.index == 0L,] <- 0

  return(list(q.list=q.list, q.index=q.index, beta.tplus=tvals$beta.tplus,
              tplus=tvals$tplus,
              KMAns=KMAns, i=seq_along(beta)))
}

## .calcSGL.coeff.smp <- function(beta, KM0, dF0, tau, time.points, RR, interval) {
##   com <- .calcSGLCoeffCommon(beta, KM0, tau, time.points)

##   return(mapply(FUN=.calcSGLBetaCoeffBasic,
##                 i=com$i, beta=beta, beta.tplus=com$beta.tplus,
##                 MoreArgs=list(q.index=com$q.index, KMAns=com$KMAns, dF0=dF0,
##                   RR=RR, interval=interval),
##                 SIMPLIFY=FALSE, USE.NAMES=FALSE))
## }

.calcSGL.coeff <- function(beta, KMFn0, KM0, dF0, p0, t0, n0, N0, n1, N1, tau,
                           time.points, RR, interval, doAnalytic) {

  com <- .calcSGLCoeffCommon(beta, KM0, tau, time.points)


  coeffs <- mapply(FUN=.calcSGLBetaCoeff,
                   i=com$i, beta=beta, beta.tplus=com$beta.tplus,
                   MoreArgs=list(q.list=com$q.list, q.index=com$q.index,
                     KMFn0=KMFn0, KMAns=com$KMAns, dF0=dF0, p0=p0, n0=n0, N0=N0,
                     n1=n1, N1=N1, RR=RR, len.F=length(KM0$t), tplus=com$tplus,
                     interval=interval, doAnalytic=doAnalytic),
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

.calcSGLInfCoeff <- function(beta, KMAns, RR, n0, N0, n1, N1, doAnalytic) {

  ans <- list()
  if(all.equal(beta, 0) == TRUE) {
    ans$Fas <- .calcSGLNegInfFas(KMAns, RR)
    if(doAnalytic)
      ans$Fas.var <- .calcSGLNegInfFasVar(KMAns, RR, n0, N0, n1, N1)
  } else {
    ans$Fas <- .calcSGLPosInfFas(KMAns, RR)
    if(doAnalytic)
      ans$Fas.var <- .calcSGLPosInfFasVar(KMAns, RR, n0, N0, n1, N1)
  }

  return(ans)
}

## .calcSGLInfCoeffBasic <- function(beta, KMAns, RR) {
##   if(beta < 0L) {
##     Fas <- .calcSGLNegInfFas(KMAns, RR)
##   } else {
##     Fas <- .calcSGLPosInfFas(KMAns, RR)
##   }

##   return(list(Fas=Fas))
## }

## .calcSGLBetaCoeffBasic <- function(i, beta, beta.tplus, q.index, KMAns, dF0,
##                                   RR, interval) {
##   if(is.infinite(beta)) {
##     return(c(list(i=i, alphahat=NA),
##              .calcSGLInfCoeffBasic(beta, KMAns, RR)))
##   }

##   alphahat <- .calc.alphahat(beta.y=beta.tplus, dF=dF0, C=RR, interval=interval)

##   if(all.equal(beta, 0) == TRUE)
##     return(list(i=i, alphahat=alphahat, Fas=KMAns$Fas))

##   w <- .calc.w(alpha=alphahat, beta.y=beta.tplus)

##   w.dF0 <- w*dF0
##   F0ai <- cumsum(w.dF0)/RR

##   Fas <- ifelse(q.index == 0L, 0L, F0ai[q.index])

##   return(list(i=i, alphahat=alphahat, Fas=Fas))
## }

.calcSGLBetaCoeff <- function(i, beta, beta.tplus, q.list, q.index,
                                 KMFn0, KMAns, dF0, p0, n0, N0, n1, N1, RR,
                                 tplus, len.F, interval, doAnalytic) {
  if(is.infinite(beta)) {
    return(c(list(i=i, alphahat=NA),
             .calcSGLInfCoeff(beta, KMAns, RR, n0, N0, n1, N1, doAnalytic)))
  }

  if(RR > 1L) {
    alphahat <- Inf
  } else {
    alphahat <- .calc.alphahat(beta.y=beta.tplus, dF=dF0, C=RR, interval=interval)}

  if(all.equal(beta, 0) == TRUE) {
    if(!doAnalytic)
      return(list(i=i, alphahat=alphahat, Fas=KMAns$Fas, FnAs=KMFn0))
    
    Fas0 <- KMAns$Fas
    FnAs <- KMFn0
    if(RR > 1L)
      w <- rep.int(1L, times=length(tplus))
    else
      w <- .calc.w(alpha=alphahat, beta.y=0)

    A1 <- -N1*p0*w*(1-w)
    B1 <- 0

    dg.nrow <- len.F + 2L
    dg <- matrix(0, ncol=length(q.index), nrow=dg.nrow,
                 dimnames=list(NULL, names(q.list)))
    
    dg[seq.int(from=0L, along=q.index)*dg.nrow + q.index + 2] <- 1
    dg <- as.data.frame(dg)
  } else {
    if(RR > 1L)
      w <- 1
    else
      w <- .calc.w(alpha=alphahat, beta.y=beta.tplus)

    w.dF0 <- w*dF0

    Fas <- cumsum(w.dF0)/RR
    FnAs <- stepfun(x=tplus, y=c(0, Fas), right=TRUE)
    
    Fas0 <- ifelse(q.index == 0, 0, Fas[q.index])

    if(!doAnalytic)
      return(list(i=i, alphahat=alphahat, Fas=Fas0, FnAs=FnAs))
    
    w.1mw.dF0 <- w.dF0*(1-w)
    sum.w.dF0 <- sum(w.dF0)
    sum.w.1mw.dF0 <- sum(w.1mw.dF0)
    diff.w <- diff(w)
  
    A1 <- sum(-N1 * w.1mw.dF0 * p0)
    B1 <- -N1*p0*-diff.w

    dg <- mapply(FUN=.calcSGL.g, q.index=q.index, q.seq=q.list,
                 MoreArgs=list(w, w.dF0, w.1mw.dF0,
                   sum.w.dF0, sum.w.1mw.dF0, diff.w, Fas),
                 SIMPLIFY=FALSE)
  }
  
  return(list(i=i, A1=A1, B1=B1, alphahat=alphahat, Fas=Fas0, FnAs=FnAs,
              Fas.var=NULL, dg=dg))
}

sensitivitySGL <- function(z, s, d, y, v, beta, tau, time.points,
                           selection, trigger, groupings,
                           empty.principal.stratum, followup.time,
                           ci=0.95, ci.method=c("analytic", "bootstrap"),
                           na.rm=FALSE, N.boot=100L, interval=c(-100,100),
                           oneSidedTest=FALSE, twoSidedTest=TRUE,
                           verbose=getOption("verbose"), isSlaveMode=FALSE) {

  ## z - group that subject belongs to
  ## s - subject met selection cirteria
  ## d - subject had event
  ## y - time until event ocurred
  withoutCdfs <- isSlaveMode && !missing(ci.method) && is.null(ci.method)
  withoutCi <- ((!isSlaveMode && !missing(ci.method) && ci.method == "") ||
                (isSlaveMode && !(!missing(ci.method) && !is.null(ci.method) &&
                                 'analytic' %in% ci.method)))

  doFollowupMethod <- !missing(followup.time) && !is.null(followup.time)

  if(!isSlaveMode) {
    ## Not running a boot strap mode
    ## Run error checks on variables.
    ErrMsg <- c(.CheckEmptyPrincipalStratum(empty.principal.stratum),
                .CheckSelection(selection, s, empty.principal.stratum),
                .CheckGroupings(groupings),
                .CheckTrigger(trigger, d),
                .CheckTau(tau),
                .CheckLength(z=z, s=s, y=y, v=v),
                .CheckZ(z, groupings, na.rm=na.rm),
                .CheckS(s, empty.principal.stratum, na.rm=na.rm),
                .CheckY(y, s, selection, na.rm=na.rm),
                .CheckD(d=d, s=s, selection=selection, na.rm=na.rm),
                .CheckV(v=v, followup.time=followup.time, na.rm=na.rm))

    if(length(ErrMsg) > 0L)
      stop(paste(ErrMsg, collapse="\n  "))

    s <- s == selection

    if(na.rm == TRUE) {
      if(doFollowupMethod)
        naIndex <- (is.na(s) | is.na(v) | is.na(z) |
                    (s & (is.na(d) | is.na(y))))
      else
        naIndex <- (is.na(s) | is.na(z) | (s & (is.na(d) | is.na(y))))

      if(any(naIndex)) {
        z <- z[!naIndex]
        s <- s[!naIndex]
        d <- d[!naIndex]
        y <- y[!naIndex]

        if(doFollowupMethod)
          v <- v[!naIndex]
      }
    }

    GroupReverse <- FALSE
    if(empty.principal.stratum[1L] == selection) {
      z <- ifelse(z == groupings[1L], TRUE, ifelse(z == groupings[2L], FALSE, NA))
      GroupReverse <- TRUE
    } else if(empty.principal.stratum[2L] == selection)
      z <- ifelse(z == groupings[2L], TRUE, ifelse(z == groupings[1L], FALSE, NA))
    d <- d == trigger
  } else {
    GroupReverse <- groupings
  }

  if(withoutCi) 
    ci.method <- NULL
  else if(isSlaveMode)
    ci.method <- "analytic"
  else
    ci.method <- sort(unique(match.arg(ci.method, several.ok=TRUE)))
  
  ## N  - number subjects
  ## N0 - number of subjects in group 0
  ## N1 - number of subjects in group 1
  N <-length(z)
  N1 <- sum(z)
  N0 <- N-N1

  if(doFollowupMethod) {
    s <- s & v < followup.time

    z0.s1 <- !z & s
    z1.s1 <- z & s

    temp <- with(survfit(Surv(v, s) ~ z, se.fit=FALSE),{
      who <- n.event > 0L
      data.frame(strata=rep.int(seq_along(strata), times=strata)[who],
                 surv=surv[who])
    })

    
    p0 <- 1L - tail(x=temp$surv[temp$strata == 1L], n=1L)
    p1 <- 1L - tail(x=temp$surv[temp$strata == 2L], n=1L)

    n0 <- p0*N0
    n1 <- p1*N1
  } else {    
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
  }

  if(all(z0.s1 == FALSE) || all(z1.s1 == FALSE) ||
     all(d[z0.s1] == FALSE) || all(d[z1.s1] == FALSE)) {
    if(isSlaveMode) {
      return(list(SCE = logical(0)))
    } else {
      stop("No times occured in one or more of the treatment arms")
    }
  }

  RR <- p1/p0
  VE <- 1L - RR

  ## summary survfit of length of time til event for group 0 and group 1.
  KM0 <- with(summary(survfit(Surv(y[z0.s1],d[z0.s1])~1L)),
              data.frame(t=time, Fas=1L - surv, Fas.var=std.err*std.err))
  KMFn0 <- stepfun(x=KM0$t, y=c(0, KM0$Fas), right=TRUE)


  KM1 <- with(summary(survfit(Surv(y[z1.s1],d[z1.s1])~1L)), 
              data.frame(t=time, Fas=1L - surv, Fas.var=std.err*std.err))

  FnAs1 <- stepfun(x=KM1$t, y=c(0, KM1$Fas), right=TRUE)

  len.t0 <- length(KM0$t)
  
  ## total length of phi vector
  len.total <- 2L + len.t0
  
  dF0 <- diff(c(0L,KM0$Fas,1L))

  doAnalyticCi <- ('analytic' %in% ci.method && !is.null(ci) && !withoutCi)
  doBootStrapCi <- ('bootstrap' %in% ci.method && !is.null(ci) && !withoutCi)

  betaOrig <- beta
  beta <- unique(sort(beta))

  timePointsOrig <- time.points
  time.points <- unique(sort(time.points))

  tpIndex <- match(time.points, timePointsOrig)
  betaIndex <- match(beta, betaOrig)

  coeffs0 <- .calcSGL.coeff(beta=beta, KMFn0=KMFn0, KM0=KM0, dF0=dF0,
                            p0=p0, n0=n0, N0=N0, n1=n1, N1=N1, tau=tau[1],
                            time.points=time.points, RR=RR,
                            interval=interval, doAnalytic=doAnalyticCi)

  coeff1 <- .calcSGLF1(time.points=time.points, KM1)
  
  SCE.dim <- c(length(beta), length(time.points))
  SCE.length <- prod(SCE.dim)
  SCE.dimnames <- list(beta=as.character(beta),
                       time.points=as.character(time.points))
  SCE <- array(numeric(SCE.length), dim=SCE.dim, dimnames=SCE.dimnames)
  
  if(!withoutCdfs) {
    FnAs0.length <- length(beta)

    alphahat <- numeric(FnAs0.length)
    names(alphahat) <- as.character(beta)

    FnAs0 <- funVector(FnAs0.length)
    names(FnAs0) <- as.character(beta)
  }
  
  for(coeff0 in coeffs0) {
    SCE[coeff0$i,] <- coeff0$Fas - coeff1$Fas
    if(!withoutCdfs) {
      FnAs0[coeff0$i] <- coeff0$FnAs
      alphahat[coeff0$i] <- coeff0$alphahat
    }
  }

  if(withoutCdfs) return(list(SCE=SCE))

  if(GroupReverse && !isSlaveMode) {
    FnAs1 <- FnAs0
    FnAs0 <- FnAs1
  }
  
  cdfs <- list()
  
  if(withoutCi) {
    if(isSlaveMode)
      return(list(SCE=SCE, alphahat=alphahat, Fas0=FnAs0, Fas1=FnAs1))

      return(structure(list(SCE=SCE[betaIndex,tpIndex,drop=FALSE],
                            beta=betaOrig, alphahat=alphahat[betaIndex],
                            Fas0=FnAs0[betaIndex], Fas1=FnAs1),
                       class=c("sensitivity.1d", "sensitivity"),
                       parameters=list(z0=groupings[1], z1=groupings[2],
                         selected=selection, trigger=trigger)))
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

  SCE.ci.dim <- c(SCE.dim, ci.probsLen, length(ci.method))
  SCE.ci.dimnames <- c(SCE.dimnames, list(ci.probs=as.character(ci.probs),
                                         ci.method=ci.method))
  
  SCE.var.dim <- SCE.ci.dim[names(SCE.ci.dimnames) != "ci.probs"]
  SCE.var.dimnames <- SCE.ci.dimnames[names(SCE.ci.dimnames) != "ci.probs"]
  
  SCE.var <- array(numeric(0),
                   dim=SCE.var.dim,
                   dimnames=SCE.var.dimnames)  

  SCE.ci <- array(numeric(0),
                  dim=SCE.ci.dim,
                  dimnames=SCE.ci.dimnames)

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
        SCE.var[coeff0$i,,'analytic'] <- coeff0$Fas.var
      } else {
        Gamma[2,2] <- coeff0$A1
        Gamma[b.col1,2] <- coeff0$B1

        IGamma <- solve(Gamma/N)
        vartheta <- crossprod(IGamma, Omega) %*% IGamma / N

        SCE.var[coeff0$i,,'analytic'] <- sapply(coeff0$dg, function(dg) {
          (dg %*% vartheta %*% dg) }) / N + coeff1$Fas.var
      }
    }

    SCE.ci[,,,'analytic'] <- array(rep.int(SCE, ci.probsLen) +
                                   rep.int(qnorm(ci.probs),
                                           rep.int(SCE.length, ci.probsLen)) *
                                   rep.int(sqrt(SCE.var[,,'analytic']),
                                           ci.probsLen),
                                   dim=c(SCE.dim, ci.probsLen))
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
                                   time.points=time.points,
                                   groupings=GroupReverse, interval=interval,
                                   ci.method=NULL,
                                   isSlaveMode=TRUE)$SCE)
      if(verbose) cat(".")
      return(ans)
    }

    vals <- apply(rows <- do.call(rbind, lapply(integer(N.boot), FUN=bootCalc,
                                        z.seq=z.seq, nVal=N,
                                        beta=beta,
                                        tau=tau[1],
                                        time.points=time.points,
                                        current.fun=current.fun,
                                        verbose=verbose)),
                  2,
                  FUN=function(x) return(c(var(x), quantile(x, probs=ci.probs))))
    N.bootActual <- nrow(rows)
    if(verbose) cat("\n")

    SCE.var.boot = vals[1,,drop=FALSE]
    SCE.ci.boot = t(vals[-1,,drop=FALSE])
    
    dim(SCE.var.boot) <- SCE.dim
    dim(SCE.ci.boot) <- c(SCE.dim, ci.probsLen)

    SCE.var[,,'bootstrap'] <- SCE.var.boot

    SCE.ci[,,,'bootstrap'] <- SCE.ci.boot
  }

  ans <- list(SCE=SCE[betaIndex,tpIndex,drop=FALSE],
              SCE.var=SCE.var[betaIndex, tpIndex, ,drop=FALSE],
              SCE.ci=SCE.ci[betaIndex, tpIndex,,, drop=FALSE],
              beta=betaOrig, alphahat=alphahat[betaIndex],
              Fas0=FnAs0[betaIndex], Fas1=FnAs1)

  if(doBootStrapCi) {
    attr(ans, 'N.boot') <- N.boot
    attr(ans, 'N.bootActual') <- N.bootActual
  }
  
  class(ans) <- c("sensitivity.1d", "sensitivity")

  return(ans)
}
