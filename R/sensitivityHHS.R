.PrepDataObj <- function(z, s, y, GroupReverse, withoutCdfs) {
  ## Various summary numbers
  ## N : total number of records in data
  ## N0, N1           : number of records in the first and second groups
  ## n0, n1           : number of selected records in the first and second group
  ## p0, p1           : probablity of a record being selected in the frist and
  ##                       second group
  ## z0.s1, z1.s0     : index vector of the selected records in the first and
  ##                       second groups
  ## y0, y1           : y values of selected records of the first and second
  ##                       groups
  ## RR               : Ratio of probablities
  ## VE               : Vacine Efficecy
  ## F0, F1           : distribution of the first and second groups
  ## y0.uniq, y1.uniq : the sorted unique values of y in the selected records
  ##                      of the first and second groups

  z0.s1 <- !z & s
  z1.s1 <- z & s

  y0 <- y[z0.s1]
  y1 <- y[z1.s1]

  N <- length(z)
  N0 <- sum(!z)
  n0 <- sum(z0.s1)
  p0 <- n0/N0

  N1 <- sum(z)
  n1 <- sum(z1.s1)
  p1 <- n1/N1

  RR <- min(p1/p0, 1L)

  VE <- 1L - RR

  Fn0 <- ecdf(y0)
  Fn1 <- ecdf(y1)

  
  c(list(withoutCdfs=withoutCdfs, GroupReverse=GroupReverse,
         z=z, s=s, y=y, z0.s1=z0.s1, z1.s1=z1.s1, N=N,
         N0=N0, N1=N1, n0=n0, n1=n1, p0=p0, p1=p1, RR=RR, VE=VE, Fn0=Fn0,
         Fn1=Fn1, y0=y0, y1=y1),
    eval(expression(list(y0.uniq=x, F0=y)), envir=environment(Fn0)),
    eval(expression(list(y1.uniq=x, F1=y)), envir=environment(Fn1)))
}

.CalcDeltaMu <- function(mu0, mu1, GroupReverse) {
  if(GroupReverse) {
    mu0 - mu1
  } else {
    mu1 - mu0
  }
}
  
.CreateACERetValHHS <- function(Fas0, obj) {
  dFas0 <- diff(c(0L, Fas0))

  mu0 <- sum(obj$y0.uniq * dFas0)
  mu1 <- mean(obj$y1)

  ACE <- .CalcDeltaMu(mu0=mu0, mu1=mu1, GroupReverse=obj$GroupReverse)

  if(obj$withoutCdfs)
    return(list(ACE=ACE))
  
  return(list(ACE=ACE,
              Fas0=Fas0, FnAs0=stepfun(x=obj$y0.uniq, y=c(0, Fas0))))
}

.CalcUpperACEHHS <- function(obj) {
  ## Variables
  ## qc : y0 value that represents the RRth percentila of y0
  ## Fas0, Fas1 : distribution funtion for the always selected stratum for the first and second groups
  ## dFas0, dFas1 : delta of Fas0 and Fas1

  qc <- quantile(obj$y0, probs=obj$RR)
  Fas0 <- ifelse(obj$y0.uniq <= qc & obj$F0 < obj$RR, obj$F0/obj$RR, 1L)
  
  .CreateACERetValHHS(Fas0=Fas0, obj=obj)
}
  
.CalcLowerACEHHS <- function(obj) {
  ## Variables
  ## qc : y0 value that represents the RRth percentila of y0
  ## Fas0, Fas1 : distribution funtion for the always selected stratum for the first and second groups
  ## dFas0, dFas1 : delta of Fas0 and Fas1

  qc <- quantile(obj$y0, probs=obj$VE)
  Fas0 <- ifelse(obj$y0.uniq >= qc, (obj$F0 - obj$VE)/obj$RR, 0)
  
  .CreateACERetValHHS(Fas0=Fas0, obj=obj)
}  

sensitivityHHS <- function(z, s, y, bound=c("upper","lower"),
                           selection, groupings, empty.principal.stratum,
                           ci=0.95, ci.method=c("bootstrap", "analytic"),
                           na.rm=FALSE, N.boot=100, oneSidedTest = FALSE,
                           twoSidedTest = TRUE, isSlaveMode=FALSE)
{
  withoutCdfs <- isSlaveMode && !missing(ci.method) && is.null(ci.method)
  withoutCi <- ((!isSlaveMode && !missing(ci.method) && ci.method == "") ||
                (isSlaveMode && !(!missing(ci.method) &&
                                 !is.null(ci.method) &&
                                 'analytic' %in% ci.method)))
  
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

    if(any(is.na(z) | is.na(s)))
      stop("s, z cannot contain any NA values")
    
    if(any(s & is.na(y)))
      stop("selected y values cannot be NA")

    GroupReverse <- FALSE
    if(empty.principal.stratum[1] == selection) {
      z <- z == groupings[1]
      GroupReverse <- TRUE
    } else if(empty.principal.stratum[2] == selection)
      z <- z == groupings[2]
  } else {
    GroupReverse <- groupings
  }
  
  bound.orig <- bound
  bound <- unique(match.arg(bound, several.ok=TRUE))
  boundIndex <- match(bound.orig, bound)

  if(withoutCi)
    ci.method <- NULL
  else if(isSlaveMode)
    ci.method <- "analytic"
  else if(missing(ci.method) || is.null(ci.method)) {
    ci.method <- 'bootstrap'
  } else {
    ci.method <- sort(unique(match.arg(ci.method, several.ok=TRUE)))
  }

  UpperIndex <- bound == "upper"
  LowerIndex <- bound == "lower"

  DoUpper <- any(UpperIndex)
  DoLower <- any(LowerIndex)
  
  datObj <- .PrepDataObj(z=z, s=s, y=y, GroupReverse=GroupReverse,
                         withoutCdfs=withoutCdfs)

  ACE.dim <- length(bound)
  ACE.length <- prod(ACE.dim)
  ACE.dimnames <- bound

  ACE <- numeric(ACE.dim)
  names(ACE) <- ACE.dimnames

  if(!withoutCdfs) {
    FnAs0 <- funVector(ACE.length)
    names(FnAs0) <- ACE.dimnames
    FnAs1 <- FnAs0
  }
  
  if(DoUpper) {
    UpperObj <- .CalcUpperACEHHS(datObj)

    ACE['upper'] <- UpperObj$ACE

    if(!withoutCdfs)
      FnAs0['upper'] <- UpperObj$FnAs0
  }

  if(DoLower) {
    LowerObj <- .CalcLowerACEHHS(datObj)

    ACE['lower'] <- LowerObj$ACE

    if(!withoutCdfs)
      FnAs0['lower'] <- LowerObj$FnAs0
  }

  if(withoutCdfs)
    return(list(ACE=ACE))

  if(!isSlaveMode && GroupReverse) {
    Fas0 <- FnAs1
    Fas1 <- FnAs0
  } else {
    Fas0 <- FnAs0
    Fas1 <- FnAs1
  }
  
  if(withoutCi) {
    return(list(ACE=ACE, Fas0=Fas0, Fas1=Fas1))
  }

  ACE.var.dim <- c(ACE.dim, length(ci.method))
  ACE.var.length <- prod(ACE.var.dim)
  ACE.var.dimnames <- c(list(bound=ACE.dimnames), list(ci.method=ci.method))

  ACE.var <- array(numeric(ACE.var.length), dim=ACE.var.dim,
                   dimnames=ACE.var.dimnames)

  ## Do analytic method
  if('analytic' %in% ci.method) {
    stop("Analytic method currently does not exist")
  }

  if(isSlaveMode) {
    return(list(ACE=ACE, ACE.var=ACE.var, Fas0=Fas0, Fas1=Fas1))
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

  ci.probs <- sort(ci.probs)
  ci.probsLen <- length(ci.probs)

  ACE.ci.dim <- c(ACE.dim, ci.probsLen, length(ci.method))
  ACE.ci.length <- prod(ACE.ci.dim)
  ACE.ci.dimnames <- c(list(bound=ACE.dimnames),
                       list(ci.probs= sprintf("%s%%",
                              as.character(ci.probs*100)),
                            ci.method=ci.method))

  ACE.ci <- array(numeric(ACE.ci.length), dim=ACE.ci.dim,
                  dimnames=ACE.ci.dimnames)

  ## run bootstrap method
  if(any(ci.method == 'bootstrap')) {
    bootACECalc <- function(i, z.seq, nVal, z, s, y, bound, GroupReverse,
                         current.fun) {
      index <- sample(z.seq, nVal, replace=TRUE)
      ans <- current.fun(z=z[index], s=s[index], y=y[index], bound=bound,
                         groupings=GroupReverse, ci.method=NULL,
                         isSlaveMode=TRUE)$ACE

      return(ans)
    }
    
    ACE.vals <- apply(do.call(rbind, lapply(integer(N.boot), z.seq=seq_along(z),
                                            nVal=datObj$N, z=datObj$z,
                                            s=datObj$s, y=datObj$y, bound=bound,
                                            GroupReverse=GroupReverse,
                                            current.fun=sys.function(),
                                            FUN=bootACECalc)),
                      2L,
                      FUN=function(x) c(var(x), quantile(x, probs=ci.probs)))

    ACE.var[,'bootstrap'] <- ACE.vals[1L,, drop=FALSE]
    ACE.ci[,,'bootstrap']  <- t(ACE.vals[-1L,, drop=FALSE])
  }

  
  if('analytic' %in% ci.method) {
    calculateCi <- function(i, norm, ACE, sqrt.ACE.var) {
      ACE[i] + norm * sqrt.ACE.var[i]
    }


    ACE.ci[,,'analytic'] <- outer(seq_along(ACE), qnorm(ci.probs),
                                  FUN=calculateCi, ACE=ACE,
                                  sqrt.ACE.var=sqrt(ACE.var[,'analytic']))
  }

  ans <- structure(list(ACE=ACE[boundIndex],
                        ACE.ci=ACE.ci[boundIndex,,, drop=FALSE],
                        ACE.var=ACE.var[boundIndex,, drop=FALSE],
                        Fas0=Fas0[boundIndex], Fas1=Fas1[boundIndex]),
                   class=c("sensitivity.0d", "sensitivity"),
                   parameters=list(z0=groupings[1], z1=groupings[2],
                     selected=selection, s0=empty.principal.stratum[1],
                     s1=empty.principal.stratum[2]))

  if('bootstrap' %in% ci.method)
    attr(ans, 'N.boot') <- N.boot

  return(ans)
}

