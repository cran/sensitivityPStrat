.PrepDataObj <- function(z, s, y, GroupReverse) {
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

  
  c(list(GroupReverse=GroupReverse,z=z, s=s, y=y, z0.s1=z0.s1, z1.s1=z1.s1, N=N,
         N0=N0, N1=N1, n0=n0, n1=n1, p0=p0, p1=p1, RR=RR, VE=VE, Fn0=Fn0,
         Fn1=Fn1),
    eval(expression(list(y0.uniq=x, F0=y)), envir=environment(ecdfF0)),
    eval(expression(list(y1.uniq=x, F1=y)), envir=environment(ecdfF1)))
}

.CalcDeltaMu <- function(mu0, mu1, GroupReverse) {
  if(GroupReverse) {
    mu0 - mu1
  } else {
    mu1 - mu0
  }
}
  
.CreateACERetValHHS <- function(Fas0, obj) {
  dFas0 <- diff(c(0L, Fas0, 1L))

  mu0 <- sum(obj$y0.uniq * dFas0)
  mu1 <- mean(obj$y1)

  ACE <- .CalcDeltaMu(mu0=mu0, mu1=mu1, GroupReverse=obj$GroupReverse)
  
  list(ACE=ACE,
       Fas0=Fas0, FnAs0=stepfun(x=obj$y0, y=c(0, Fas0)))
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
                           selection, groupings, empty.principle.stratum,
                           ci=0.95, ci.method=c("bootstrap", "analytic"),
                           na.rm=FALSE, N.boot=100, oneSidedTest = FALSE,
                           twoSidedTest = TRUE)
{
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
    isEmptyPrincipleStratumMissing <- missing(empty.principle.stratum)
    
    if(isEmptyPrincipleStratumMissing) 
      ErrMsg <- c(ErrMsg, "'empty.principle.stratum' argument must be supplied")
    else {
      if(is.null(empty.principle.stratum) ||
         length(empty.principle.stratum) != 2L)
      ErrMsg <- c(ErrMsg,
                  "'empty.principle.stratum' argument must be a two element vector")

      if(length(empty.principle.stratum) &&
         any(is.na(empty.principle.stratum)))
      ErrMsg <- c(ErrMsg,
                  "'empty.principle.stratum' may not contain a NA")

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
      
      if(length(selection) > 0L && !isEmptyPrincipleStratumMissing &&
         !all((selection %in% empty.principle.stratum) | is.na(selection)))
        ErrMsg <- c(ErrMsg,
                    "'selection' value does not appear in specified levels of 's' from 'empty.principle.stratum'")
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
      
      if(!isEmptyPrincipleStratumMissing &&
         !all(s %in% empty.principle.stratum | is.na(s)))
        ErrMsg <- c(ErrMsg,
                    "All values of 's' must match one of the two values in 'empty.principle.stratum'")

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

    ParamConflict <- c(!missing(psi), !missing(Pi), !missing(phi))
    if(sum(ParamConflict) > 1)
      ErrMsg <- c(ErrMsg,
                  "only one of 'psi', 'Pi', or 'phi' can be specified")
    else if(sum(ParamConflict < 1))
      ErrMsg <- c(ErrMsg,
                  "one of 'psi', 'Pi', or 'phi' must be specified")
      
    
    if(length(ErrMsg) > 0L)
      stop(paste(ErrMsg, collapse="\n  "))

    s <- s == selection

    GroupReverse <- FALSE
    if(empty.principle.stratum[1] == selection) {
      z <- z == groupings[1]
      GroupReverse <- TRUE
    } else if(empty.principle.stratum[2] == selection)
      z <- z == groupings[2]
  } else {
    GroupReverse <- ci
  }
  
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

  bound <- sort(unique(match.arg(bound, several.ok=TRUE)))
  
  ci.method <- match.arg(ci.method)

  UpperIndex <- bound == "upper"
  LowerIndex <- bound == "lower"

  DoUpper <- any(UpperIndex)
  DoLower <- any(LowerIndex)
  
  datObj <- .PrepDataObj(z=z, s=s, y=y, GroupReverse=GroupReverse)

  
  if(doUpper)
    UpperObj <- .CalcUpperACEHHS(obj)

  if(doLower)
    LowerObj <- .CalcLowerACEHHS(obj)

  if(doUpper & doLower) {
    ACE <- c(lower=LowerObj$ACE, upper=UpperObj$ACE)
    FnAs0 <- list(lower=LowerObj$FunAs0, upper=UpperObj$FunAs0)
  } else if(doUpper) {
    ACE <- c(upper=UpperObj$ACE)
    FnAs0 <- list(upper=UpperObj$FunAs0)
  } else if(doLower) {
    ACE <- c(lower=LowerObj$ACE)
    FnAs0 <- list(lower=LowerObj$FunAs0)
  }

  if(isSlaveMode) {
    return(list(ACE=ACE))
  }

  if(GroupReverse) {
    cdfs <- list(Fas0=list(obj$Fn1), Fas1=FnAs0)
  } else {
    cdfs <- list(Fas0=FnAs0, Fas1=list(obj$Fn1))
  }

  if(is.null(ci)) {
    return(list(ACE=ACE), cdfs)
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

  ## run bootstrap method
  ACE.info <- lapply(ci.method, switch,
                     analytic = stop("Analytic method currently does not exist"),
                     bootstrap = {
                       ACE.list <- matrix(nrow=N.boot, ncol=length(ACE),
                                          dimnames=list(NULL, bound.uniq))
                       index <- seq_len(N)
                       for(i in seq_len(N.boot)) {
                         new.index <- sample(index, N, replace=TRUE)
                         ACE.list[i,] <- Recall(z=z[new.index],s=s[new.index],
                                                y=y[new.index],
                                                bound=bound.uniq,
                                                ci=GroupReverse,
                                                ci.method=NULL,
                                                selection=TRUE,
                                                groupings=c(FALSE,TRUE),
                                                empty.principle.stratum=c(FALSE,TRUE))$ACE
                       }
                       list(ACE.var=diag(var(ACE.list))[bIndex],
                            ACE.ci=t(apply(ACE.list, MARGIN=2, FUN=quantile,
                              probs=ci.probs))[bIndex,])
                     }
                     )

  names(ACE.info) <- ci.method

  DIM <- function(x) if(is.array(x) || is.data.frame(x)) dim(x) else length(x)
  DIMNAMES <- function(x) if(is.array(x) || is.data.frame(x)) dimnames(x) else list(names(x))
  bind <- function(name, ...) {
    dotlist <- list(...)
    upper <- length(dotlist)
    uppernames <- list(names(dotlist))

    dims <- DIM(dotlist[[1]])
    dnames <- DIMNAMES(dotlist[[1]])

    array(c(...), dim=c(dims, upper), dimnames=c(dnames, uppernames))
  }
  ACE.info <- do.call(function(...) mapply(FUN=bind, ..., USE.NAMES=TRUE,
                                           SIMPLIFY=FALSE),
                      ACE.info)

  ans <- structure(c(list(ACE=ACE[bIndex], bound=bound[bIndex]),
                     cdfs, ACE.info),
                   class=c("sensitivity2d", "sensitivity"),
                   parameters=list(z0=groupings[1], z1=groupings[2],
                     selected=selection, s0=empty.principle.stratum[1],
                     s1=empty.principle.stratum[2]))
  if('bootstrap' %in% ci.method)
    attr(ans, 'N.boot') <- N.boot

  return(ans)
}

