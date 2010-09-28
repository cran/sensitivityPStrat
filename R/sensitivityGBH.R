sensitivityGBH <- function(z, s, y, beta, selection, groupings,
                           empty.principal.stratum, ci=0.95,
                           ci.method=c("analytic", "bootstrap"), na.rm=FALSE, 
                           N.boot=100, oneSidedTest = FALSE,
                           twoSidedTest = TRUE)
{

  doInfinite <- any(is.infinite(beta))
  doFinite <- any(is.finite(beta))

  if(doInfinite && !doFinite) {
    funCall <- match.call()
    funCall$beta <- c("upper", "lower")[match(beta, c(-Inf, Inf), nomatch=0)]
    names(funCall)[names(funCall) == "beta"] <- "bound"

    return(eval(c))
  }

  if(na.rm == TRUE) {
    naIndex <- !(is.na(s) | is.na(z))

    z <- z[naIndex]
    s <- s[naIndex]
    y <- y[naIndex]
  }
  

  if(!missing(ci.method) && is.null(ci.method))
    isSlaveMode <- TRUE
  else
    isSlaveMode <- FALSE

  if(!isSlaveMode) {
    ## Not running a boot strap mode
    ## Run error checks on variables.
    ErrMsg <- character(0)

    isZMissing <- missing(z)
    isSMissing <- missing(s)
    isYMissing <- missing(y)

    isSelectionMissing <- missing(selection)
    isGroupingsMissing <- missing(groupings)
    isEmptyPrincipalStratumMissing <- missing(empty.principal.stratum)
    
    if(isEmptyPrincipalStratumMissing) 
      ErrMsg <- c(ErrMsg, "'empty.principal.stratum' argument must be supplied")
    else {
      if(is.null(empty.principal.stratum) ||
         length(empty.principal.stratum) != 2L)
      ErrMsg <- c(ErrMsg,
                  "'empty.principal.stratum' argument must be a two element vector")

      if(length(empty.principal.stratum) &&
         any(is.na(empty.principal.stratum)))
      ErrMsg <- c(ErrMsg,
                  "'empty.principal.stratum' may not contain a NA")

    }

    if(isSelectionMissing)
      ErrMsg <- c(ErrMsg, "'selection' argument must be supplied")
    else {
      if(length(selection) > 0L && any(is.na(selection)))
        ErrMsg <- c(ErrMsg,
                    "'selection' may not be NA")

      if(length(selection) != 1L)
        ErrMsg <- c(ErrMsg,
                    "'selection' argument must be a single element vector")
      
      if(length(selection) > 0L && !isEmptyPrincipalStratumMissing &&
         !all((selection %in% empty.principal.stratum) | is.na(selection)))
        ErrMsg <- c(ErrMsg,
                    "'selection' value does not appear in specified levels of 's' from 'empty.principal.stratum'")
    }

    if(isGroupingsMissing)
      ErrMsg <- c(ErrMsg,
                  "'groupings' argument must be supplied")
    else {
      if(is.null(groupings) || length(groupings) != 2L)
        ErrMsg <- c(ErrMsg,
                    "'groupings' argument must be a two element vector")

      if(length(groupings) && any(is.na(groupings)))
        ErrMsg <- c(ErrMsg,
                    "'groupings' may not contain a NA")
    }
    
    vectorLength <- NULL
    SameLength <- TRUE

    if(isZMissing)
      ErrMsg <- c(ErrMsg,
                  "'z' argument must be supplied")
    else {
      if(any(is.na(z)) && na.rm)
        ErrMsg <- c(ErrMsg,
                    "'z' cannot contain any NA values")

      if(!isGroupingsMissing && !all(z %in% groupings | is.na(z)))
        ErrMsg <- c(ErrMsg,
                    "All values of 'z' must match one of the two values in 'groupings'")

      if(is.null(vectorLength))
        vectorLength <- length(z)
    }

    if(isSMissing)
      ErrMsg <- c(ErrMsg,
                  "'s' argument must be supplied")
    else {
      if(any(is.na(s)) && !na.rm)
        ErrMsg <- c(ErrMsg,
                    "argument 's' cannot contain any NA values")
      
      if(!isEmptyPrincipalStratumMissing &&
         !all(s %in% empty.principal.stratum | is.na(s)))
        ErrMsg <- c(ErrMsg,
                    "All values of 's' must match one of the two values in 'empty.principal.stratum'")

      if(is.null(vectorLength))
        vectorLength <- length(s)
      else if(SameLength == TRUE && vectorLength != length(s))
        SameLength <- FALSE
    }

    if(isYMissing) 
      ErrMsg <- c(ErrMsg,
                  "'y' argument must be supplied")
    else {
      if(is.null(vectorLength))
        vectorLength <- length(y)
      else if(SameLength == TRUE && vectorLength != length(y))
        SameLength <- FALSE

      if(!isSelectionMissing && !isSMissing && length(s) == length(y) &&
         any(s %in% selection && !is.na(s) && is.na(y)))
         ErrMsg <- c(ErrMsg,
                     sprintf("argument 'y' cannont contain a NA value if the corresponding 's' is %s",
                          paste(selection, collapse=",")))
    }

    if(SameLength == FALSE)
        ErrMsg <- c(ErrMsg,
                    "'z', 's', 'y' are not all the same length")
    
    if(length(ErrMsg) > 0L)
      stop(paste(ErrMsg, collapse="\n  "))

    s <- s == selection

    GroupReverse <- FALSE
    if(empty.principal.stratum[1] == selection) {
      z <- z == groupings[1]
      GroupReverse <- TRUE
    } else if(empty.principal.stratum[2] == selection)
      z <- z == groupings[2]
  } else {
    GroupReverse <- ci
  }

  ci.method <- sort(unique(match.arg(ci.method, several.ok=TRUE)))
  n.method <- length(ci.method)

  if(na.rm == TRUE) {
    naIndex <- !(is.na(s) | is.na(z) | (s & is.na(y)))
    z <- z[naIndex]
    s <- s[naIndex]
    y <- y[naIndex]
  }

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

  params <- .calc.ecdf(y0)
  F0 <- params$F
  dF0 <- diff(c(0, F0))
  y0.uniq <- params$vals
  
  params <- .calc.ecdf(y[z1.s1])
  Fas1 <- params$F
  dFas1 <- diff(c(0, Fas1))
  Fas1 <- matrix(Fas1, ncol=1)
  y1.uniq <- params$vals

  mu1 <- mean(y[z1.s1], na.rm=TRUE)

  w.calc <- function(alpha, beta, y)
    1/(1 + exp(-alpha - beta*y))

  alpha.est <- function(alpha, beta, y, dF, C)
    (sum(w.calc(alpha=alpha, beta=beta, y=y) * dF) - C)^2

  calc.dFas0AndAlpha <- function(beta, y, dF, C) {
    alphahat <- optimize(f=alpha.est, interval=c(-100,100), beta=beta,
                         y=y, dF=dF, C=C)$minimum

    if(alphahat > 90 || alphahat < -90) {
      warning("optimize overflow alphahat value invalid")
    }
    
    w <- w.calc(alpha=alphahat, beta=beta, y=y)
    
    return(c(alphahat=alphahat, dF*w/C))
  }

  Fas0 <- matrix(NA, nrow=length(y0.uniq), ncol=length(beta))
  alphahat <- numeric(length(beta))
  ACE.dim <- length(beta)
  ACE.dimnames <- format(beta, trim=TRUE)

  ACE <- numeric(ACE.dim)
  names(ACE) <- ACE.dimnames
  
  if(doFinite) {
    mu0 <- numeric(length(finiteBeta))

    temp <- sapply(finiteBeta, FUN=calc.dFas0AndAlpha, y=y0.uniq, dF=dF0, C=RR)
    
    alphahat[finiteIndex] <- temp[1,]
    dFas0 <- temp[-1,, drop=FALSE]

    mu0 <- colSums(y0.uniq * dFas0)

    Fas0[finiteIndex] <- lapply(X=as.data.frame(dFas0), FUN=function(x) {
      x <- cumsum(x)
      stepfun(y0.uniq, y=c(0, x))
    })

    if(GroupReverse)
      ACE[finiteIndex] <- mu0 - mu1
    else
      ACE[finiteIndex] <- mu1 - mu0
  }

  if(doInfinite) {
    if(length(ci.method[ci.method != "bootstrap"]) == 0) {
      infiniteACE.info <- sensitivityHHS(z=z,s=s,y=y, bound=bounds,
                                         selection=TRUE, groupings=c(FALSE, TRUE),
                                         empty.principal.stratum=c(FALSE, TRUE),
                                         ci=GroupReverse, ci.method=NULL)
    } else {
      infiniteACE.info <- sensitivityHHS(z=z,s=s,y=y, bound=bounds,
                                         selection=TRUE, groupings=c(FALSE, TRUE),
                                         empty.principal.stratum=c(FALSE,TRUE),
                                         ci=ci, ci.method='analytic')
    }

    if(GroupReverse) 
      Fas0[infiniteIndex] <- infiniteACE.info$Fas1
    else
      Fas0[infiniteIndex] <- infiniteACE.info$Fas0
    
    ACE[infiniteIndex] <- infiniteACE.info$ACE

    alphahat[infiniteIndex] <- NA
  }

  cdfs <- list(beta=beta[bIndex], alphahat=alphahat[bIndex])
  if(GroupReverse)
    cdfs <- c(cdfs, list(y0=y1.uniq, Fas0=Fas1, y1=y0.uniq, Fas1=Fas0[bIndex]))
  else
    cdfs <- c(cdfs, list(y0=y0.uniq, Fas0=Fas0[bIndex], y1=y1.uniq, Fas1=Fas1))

  if(is.null(ci) || isSlaveMode) {
    return(c(list(ACE=ACE), cdfs))
  }
  
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

  ACE.ci.dim <- c(ACE.dim, length(ci.method), ci.probsLen)
  ACE.ci.length <- prod(ACE.ci.dim)
  ACE.ci.dimnames <- c(list(ACE.dimnames),
                       list(ci.method=ci.method,
                            ci.probs= sprintf("%s%%",
                              as.character(ci.probs*100))))

  ACE.ci <- array(numeric(ACE.ci.length), dim=ACE.ci.dim,
                  dimnames=ACE.ci.dimnames)

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
      ACE.ci[infiniteIndex,'analytic',] <- infiniteACE.info$ACE.ci[,'analytic',]
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

    calculateCi <- function(i, norm, ACE, sqrt.ACE.var) {
      ACE[i] + norm * sqrt.ACE.var[i]
    }

    ACE.ci[,'analytic',] <- outer(seq_along(ACE), qnorm(ci.probs),
                                  FUN=calculateCi, ACE=ACE,
                                  sqrt.ACE.var=sqrt(ACE.var[,'analytic']))
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
                                                  bound=bounds, selection=TRUE,
                                                  groupings=c(FALSE, TRUE),
                                                  empty.principal.stratum=c(FALSE, TRUE),
                                                  ci=GroupReverse,
                                                  ci.method=NULL)$ACE
      
      ACE.vals[finiteIndex] <- current.fun(z=new.z, s=new.s, y=new.y,
                                           beta=beta[finiteIndex],
                                           selection=TRUE,
                                           groupings=c(FALSE, TRUE),
                                           empty.principal.stratum=c(FALSE, TRUE),
                                           ci=GroupReverse,
                                           ci.method=NULL)$ACE

      ACE.vals
    }, simplify=FALSE), recursive=FALSE)

    ans <- sapply(split(ACE.list, rep.int(beta, times=N.boot)),
                  function(ACE.vals, probs) c(var(ACE.vals),
                                              quantile(ACE.vals, probs=probs)),
                  probs=ci.probs)
    ACE.ci[,'bootstrap',] <- t(ans[-1,])
    ACE.var[,'bootstrap'] <- ans[1,]
  }


  ans <- structure(c(list(ACE=ACE[bIndex],
                          ACE.ci=ACE.ci[bIndex,,,drop=FALSE],
                          ACE.var=ACE.var[bIndex,, drop=FALSE]), cdfs),
                   class=c("sensitivity2d", "sensitivity"),
                   parameters=list(z0=groupings[1], z1=groupings[2],
                     selected=selection, s0=empty.principal.stratum[1],
                     s1=empty.principal.stratum[2]))
  if('bootstrap' %in% ci.method)
    attr(ans, 'N.boot') <- N.boot
  
  return(ans)
}
