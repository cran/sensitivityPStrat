sensitivityGBH <- function(z, s, y, beta, selection, groupings,
                           empty.principal.stratum, ci=0.95,
                           ci.method=c("analytic", "bootstrap"), na.rm=FALSE, 
                           N.boot=100, interval=c(-100,100),
                           oneSidedTest = FALSE, twoSidedTest = TRUE,
                           isSlaveMode=FALSE)
{
  withoutCdfs <- isSlaveMode && !missing(ci.method) && is.null(ci.method)
  withoutCi <- ((!isSlaveMode && !missing(ci.method) && ci.method == "") ||
                (isSlaveMode && !(!missing(ci.method) &&
                                 !is.null(ci.method) &&
                                 'analytic' %in% ci.method)))

  doInfinite <- any(is.infinite(beta))
  doFinite <- any(is.finite(beta))

  ## if(doInfinite && !doFinite) {
  ##   funCall <- match.call()
  ##   funCall$beta <- c("upper", "lower")[match(beta, c(-Inf, Inf), nomatch=0)]
  ##   names(funCall)[names(funCall) == "beta"] <- "bound"
  ##   funCall[[1]] <- as.name("sensitivityHHS")

  ##   return(eval(funCall, envir=parent.frame()))
  ## }

  if(!isSlaveMode) {
    ## Not running a boot strap mode
    ## Run error checks on variables.
    ErrMsg <- c(.CheckEmptyPrincipalStratum(empty.principal.stratum),
                .CheckSelection(selection, s, empty.principal.stratum),
                .CheckGroupings(groupings),
                .CheckLength(z=z, s=s, y=y),
                .CheckZ(z, groupings, na.rm=na.rm),
                .CheckS(s, empty.principal.stratum, na.rm=na.rm),
                .CheckY(y, s, selection, na.rm=na.rm))

    if(length(ErrMsg) > 0L)
      stop(paste(ErrMsg, collapse="\n  "))

    s <- s == selection

    if(na.rm == TRUE) {
      naIndex <- !(is.na(s) | is.na(z) | (s & is.na(y)))

      z <- z[naIndex]
      s <- s[naIndex]
      y <- y[naIndex]
    }
    
    GroupReverse <- FALSE
    if(empty.principal.stratum[1] == selection) {
      z <- z == groupings[1]
      GroupReverse <- TRUE
    } else if(empty.principal.stratum[2] == selection)
      z <- z == groupings[2]
  } else {
    GroupReverse <- groupings
  }

  if(withoutCi)
    ci.method <- NULL
  else if(isSlaveMode)
    ci.method <- "analytic"
  else
    ci.method <- sort(unique(match.arg(ci.method, several.ok=TRUE)))

  n.method <- length(ci.method)

  if(any(is.na(z) | is.na(s)))
    stop("s, z cannot contain any NA values")
  
  if(any(s & is.na(y)))
    stop("selected y values cannot be NA")
  
  beta.orig <- beta
  beta <- unique(beta)
  bIndex <- match(beta.orig, beta)

  finiteBeta <- beta[is.finite(beta)]
  finiteIndex <- match(finiteBeta, beta, nomatch=0)
  infiniteBeta <- beta[is.infinite(beta)]
  infiniteIndex <- match(infiniteBeta, beta, nomatch=0)

  bounds <- c("upper", "lower")[match(infiniteBeta, c(-Inf, Inf), nomatch=0)]

  z0.s1 <- !z & s
  z1.s1 <- z & s
  
  y0 <- y[z0.s1]
  
  N <- length(z)
  N0 <- sum(!z)
  n0 <- sum(z0.s1)
  p0 <- n0/N0

  N1 <- sum(z)
  n1 <- sum(z1.s1)
  p1 <- n1/N1

  RR <- min(p1/p0, 1)

  Fas0 <- ecdf(y0)
  y0.uniq <- knots(Fas0)
  dF0 <- diff(c(0, Fas0(y0.uniq)))
  
  Fas1 <- ecdf(y[z1.s1])
  y1.uniq <- knots(Fas1)
  dFas1 <- diff(c(0, Fas1(y1.uniq)))
  
  mu1 <- mean(y[z1.s1], na.rm=TRUE)

  calc.dFas0AndAlpha <- function(beta, y, dF, C, interval) {
    alphahat <- .calc.alphahat(beta.y=beta*y, dF=dF, C=C, interval=interval)
    
    w <- .calc.w(alpha=alphahat, beta=beta*y)
    
    return(c(alphahat=alphahat, dF*w/C))
  }

  if(!withoutCdfs)
    Fas0 <- funVector(length(beta))
  
  alphahat <- numeric(length(beta))
  ACE.dim <- length(beta)
  ACE.dimnames <- format(beta, trim=TRUE)

  ACE <- numeric(ACE.dim)
  names(ACE) <- ACE.dimnames
  
  if(doFinite) {
    mu0 <- numeric(length(finiteBeta))

    temp <- sapply(finiteBeta, FUN=calc.dFas0AndAlpha, y=y0.uniq, dF=dF0, C=RR,
                   interval=interval)
    
    alphahat[finiteIndex] <- temp[1,]
    dFas0 <- temp[-1,, drop=FALSE]

    mu0 <- colSums(y0.uniq * dFas0)

    if(!withoutCdfs) {
      Fas0[finiteIndex] <- lapply(X=as.data.frame(dFas0), FUN=function(x) {
        x <- cumsum(x)
        stepfun(y0.uniq, y=c(0, x))
      })
    }

    if(GroupReverse)
      ACE[finiteIndex] <- mu0 - mu1
    else
      ACE[finiteIndex] <- mu1 - mu0
  }

  if(doInfinite) {
    if(withoutCdfs) {
      infiniteACE.info <- sensitivityHHS(z=z,s=s,y=y, bound=bounds,
                                         groupings=GroupReverse,
                                         ci.method=NULL, isSlaveMode=TRUE)
    } else {
      infiniteACE.info <- sensitivityHHS(z=z,s=s,y=y, bound=bounds,
                                         groupings=GroupReverse,
                                         ci.method=ci.method, isSlaveMode=TRUE)

      if(!withoutCdfs) {
        Fas0[infiniteIndex] <- infiniteACE.info$Fas0
      }
      
      alphahat[infiniteIndex] <- NA
    }

    ACE[infiniteIndex] <- infiniteACE.info$ACE
  }

  if(withoutCdfs) {
    return(list(ACE=ACE))
  }
  
  if(withoutCi) {
    if(!isSlaveMode && GroupReverse)
      cdfs <- list(Fas0=Fas1, Fas1=Fas0, alphahat=alphahat)
    else
      cdfs <- list(Fas0=Fas0, Fas1=Fas1, alphahat=alphahat)
    
    return(c(list(ACE=ACE), cdfs))
  }

  if(!isSlaveMode) {
    if(twoSidedTest) {    
      ci.probs <- unlist(lapply(ci, FUN=function(ci) {
        if(ci < 0.5)
          c(ci, 1L) - ci/2L
        else
          c(0L, ci) + (1-ci)/2
      }))
    } else {
      ci.probs <- NULL
    }

    if(oneSidedTest) {
      ci.probs <- c(ci.probs, ci)
    }

    ci.probs <- unique(sort(ci.probs))
    ci.probsLen <- length(ci.probs)

    ACE.ci.dim <- c(ACE.dim, ci.probsLen, length(ci.method))
    ACE.ci.length <- prod(ACE.ci.dim)
    ACE.ci.dimnames <- c(list(ACE.dimnames),
                         list(ci.probs= sprintf("%s%%",
                                as.character(ci.probs*100)),
                              ci.method=ci.method))

    ACE.ci <- array(numeric(ACE.ci.length), dim=ACE.ci.dim,
                    dimnames=ACE.ci.dimnames)
  }
  
  ACE.var.dim <- c(ACE.dim, length(ci.method))
  ACE.var.length <- prod(ACE.var.dim)
  ACE.var.dimnames <- c(list(ACE.dimnames), list(ci.method=ci.method))

  ACE.var <- array(numeric(ACE.var.length), dim=ACE.var.dim,
                   dimnames=ACE.var.dimnames)

  ## run bootstrap method
  if(any(ci.method == 'analytic')) {
    Omega <- matrix(nrow=5,ncol=5)
    names(ACE.var) <- beta
             
    if(doInfinite) {
      ACE.ci[infiniteIndex,,'analytic'] <- infiniteACE.info$ACE.ci[,,'analytic']
      ACE.var[infiniteIndex,'analytic'] <- infiniteACE.info$ACE.var[,'analytic']
    }

    for(i in finiteIndex) {
      U <- rbind((p0-s)*(1-z),
                 (p1-s)*z,
                 s*(1/(exp(-beta[i]*y-alphahat[i])+1)-p1/p0)*(1-z),
                 s*(mu0[i]-p0*y/(p1*(exp(-beta[i]*y-alphahat[i])+1)))*(1-z),
                 s*(mu1-y)*z)
      
      .sumCrossUpperTri(Omega) <- U

      Omega <- Omega/N
      Omega <- .foldUpperTri(Omega)
      
      Gamma <- matrix(colSums(cbind(1-z,
                                    0,
                                    p1*s*(1-z)/p0^2,
                                    -s*y*(1-z)/(p1*(exp(-beta[i]*y-alphahat[i])+1)),
                                    0,

                                    
                                    0,
                                    z,
                                    -s*(1-z)/p0,
                                    p0*s*y*(1-z)/(p1^2*(exp(-beta[i]*y-alphahat[i])+1)),
                                    0,

                                    
                                    0,
                                    0,
                                    s*exp(-beta[i]*y-alphahat[i])*(1-z)/(exp(-beta[i]*y-alphahat[i])+1)^2,
                                    -p0*s*y*exp(-beta[i]*y-alphahat[i])*(1-z)/(p1*(exp(-beta[i]*y-alphahat[i])+1)^2),
                                    0,

                                    
                                    0,
                                    0,
                                    0,
                                    s*(1-z),
                                    0,

                                    
                                    0,
                                    0,
                                    0,
                                    0,
                                    s*z),
                              na.rm=TRUE),
                      nrow=5, byrow=FALSE)/N
      
      IGamma <- solve(Gamma)
      vartheta <- IGamma %*% Omega %*% t(IGamma) / N
      ACE.var[i,'analytic'] <- vartheta[4,4]+vartheta[5,5] - 2*vartheta[4,5]
    }

    if(!isSlaveMode) {
      calculateCi <- function(i, norm, ACE, sqrt.ACE.var) {
        ACE[i] + norm * sqrt.ACE.var[i]
      }

      ACE.ci[,,'analytic'] <- outer(seq_along(ACE), qnorm(ci.probs),
                                    FUN=calculateCi, ACE=ACE,
                                    sqrt.ACE.var=sqrt(ACE.var[,'analytic']))
    }
  }

  if(isSlaveMode) {
    cdfs <- list(Fas0=Fas0, Fas1=Fas1, alphahat=alphahat)

    return(c(list(ACE=ACE, ACE.var=ACE.var), cdfs))
  }
  
  if(any(ci.method == 'bootstrap')) {
    index <- seq_len(N)
    current.fun <- sys.function()
    ACE.list <- unlist(replicate(N.boot, {
      ACE.vals <- numeric(length(beta))
      names(ACE.vals) <- beta
      new.index <- sample(index, N, replace=TRUE)
      new.z <- z[new.index]
      new.s <- s[new.index]
      new.y <- y[new.index]
      if(doInfinite)
        ACE.vals[infiniteIndex] <- sensitivityHHS(z=new.z, s=new.s, y=new.y,
                                                  bound=bounds,
                                                  groupings=GroupReverse,
                                                  ci.method=NULL,
                                                  isSlaveMode=TRUE)$ACE
      
      ACE.vals[finiteIndex] <- current.fun(z=new.z, s=new.s, y=new.y,
                                           beta=beta[finiteIndex],
                                           groupings=GroupReverse,
                                           ci.method=NULL, isSlaveMode=TRUE)$ACE

      ACE.vals
    }, simplify=FALSE), recursive=FALSE)

    ans <- sapply(split(ACE.list, rep.int(beta, times=N.boot)),
                  function(ACE.vals, probs) c(var(ACE.vals),
                                              quantile(ACE.vals, probs=probs)),
                  probs=ci.probs)
    ACE.ci[,,'bootstrap'] <- t(ans[-1,])
    ACE.var[,'bootstrap'] <- ans[1,]
  }

  if(!isSlaveMode && GroupReverse)
    cdfs <- list(Fas0=Fas1, Fas1=Fas0[bIndex])
  else
    cdfs <- list(Fas0=Fas0[bIndex], Fas1=Fas1)

  ans <- structure(c(list(ACE=ACE[bIndex],
                          ACE.ci=ACE.ci[bIndex,,,drop=FALSE],
                          ACE.var=ACE.var[bIndex,, drop=FALSE],
                          beta=beta[bIndex], alphahat=alphahat[bIndex]),
                     cdfs),
                   class=c("sensitivity.1.0d", "sensitivity.0d", "sensitivity"),
                   parameters=list(z0=groupings[1], z1=groupings[2],
                     selected=selection, s0=empty.principal.stratum[1],
                     s1=empty.principal.stratum[2]))
  if('bootstrap' %in% ci.method)
    attr(ans, 'N.boot') <- N.boot
  
  return(ans)
}
