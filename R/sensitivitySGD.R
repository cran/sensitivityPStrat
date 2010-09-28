.calc.coeff <- function(Pi, p0, p1, beta0, beta1, dF0, dF1, 
                       t0, t1, tau, time.points) {
  calc.time.seq <- function(time.point, times) which(times <= time.point)
  
  q.seq.list <- lapply(time.points, t0, FUN=calc.time.seq)
  r.seq.list <- lapply(time.points, t1, FUN=calc.time.seq)
  names(q.seq.list) <- names(r.seq.list) <- time.points

  tplus0 <- c(pmin(t0, tau), tau)
  tplus1 <- c(pmin(t1, tau), tau)

  calc.beta.tplus <- function(beta,i,time) return(list(bt=beta*time, i=i))
  beta.tplus0 <- mapply(FUN=calc.beta.tplus,
                        beta=beta0, i=seq_along(beta0),
                        MoreArgs=list(time=tplus0),
                        USE.NAMES=TRUE, SIMPLIFY=FALSE)
  beta.tplus1 <- mapply(FUN=calc.beta.tplus,
                        beta=beta1, i=seq_along(beta1),
                        MoreArgs=list(time=tplus1),
                        USE.NAMES=TRUE, SIMPLIFY=FALSE)

  Pi..p0 <- Pi/p0
  Pi..p1 <- Pi/p1

  coeffs <- mapply(FUN=.calc.Pi.coeff,
                   Pi..p0=Pi..p0, Pi..p1=Pi..p1, i=seq.int(along.with=Pi),
                   MoreArgs=list(beta.tplus0=beta.tplus0,
                     beta.tplus1=beta.tplus1,
                     q.seq.list=q.seq.list, r.seq.list=r.seq.list,
                     dF0=dF0, dF1=dF1),
                   USE.NAMES = FALSE, SIMPLIFY=FALSE)
  return(coeffs)
}

.calc.Pi.coeff <- function(Pi..p0, Pi..p1, i, beta.tplus0, beta.tplus1,
                          q.seq.list, r.seq.list, dF0, dF1, dV, dW, VmF0, WmF1,
                          p0, p1) {
  beta0.coeff <- lapply(beta.tplus0, q.list=q.seq.list, dF=dF0, Pi..p=Pi..p0,
                        FUN=.calc.beta.coeff)
  beta1.coeff <- lapply(beta.tplus1, q.list=r.seq.list, dF=dF1, Pi..p=Pi..p1,
                        FUN=.calc.beta.coeff)

  return(list(beta0.coeff=beta0.coeff, beta1.coeff=beta1.coeff, i=i))
}    

.calc.beta.coeff <- function(beta.tplus, q.list, dV, dF, p, Pi..p, Pi.N) {
  calc.alphahat <- function(beta.tplus, dF, B) {
    ee <- function(a, beta.tplus, dF, B) {
      tmp <- sum((1L + exp(-a - beta.tplus))^(-1) * dF) - B
      tmp*tmp
    }

    overflow <- double(0)
    range <- c(-100L,100L)
    while(TRUE) {
      alphahat <- optimize(ee, range, beta.tplus=beta.tplus, dF=dF,
                           B=B)$minimum

      limits <- range + c(-10L, 10L)
      if(alphahat < limits[2] && alphahat > limits[1])
        break

      
      warning("optimize overflow alphahat value invalid ", alphahat,
              imediate=TRUE)
      overflow <- c(overflow, alphahat)
      if(alphahat > limits[2])
        range <- range + 75
      else if(alphahat < limits[1])
        range <- range - 75
      else
        stop("limit failure")      
    }

    nInvalid <- length(overflow)
    if(nInvalid > 0L)
      warning("optimize overflow detected\n    Final alphahat value is ",
              alphahat, "\n    ", nInvalid, "invalid alphahat values found:\n",
              paste(alphahat, collapse=', '), sep='')

    return(alphahat)
  }

  alphahat <- calc.alphahat(beta.tplus$bt, dF=dF, B=Pi..p)

  w <- (1L + exp(-alphahat - beta.tplus$bt))^(-1)

  SCE <- sapply(q.list, w.dF=w*dF,
                FUN=function(q.seq, w.dF) sum(w.dF[q.seq])) / Pi..p
  
  return(list(i=beta.tplus$i, alphahat=alphahat, SCE=SCE))
}

.makeBootstrapEvntIndx <- function(s, indx.seq, N) {
  numEvents <- 0L
  index <- integer(0L)
  storage.mode(indx.seq) <- "integer"
  storage.mode(N) <- "integer"
  
  repeat {
    subIndx <- sample(indx.seq, 1000L, replace=TRUE)
    numSubEvents <- sum(s[subIndx])
    newNumEvents <- numEvents + numSubEvents

    if(newNumEvents > N)
      subIndx <- subIndx[cumsum(s[subIndx]) + numEvents <= N]

    index <- c(index, subIndx)

    if(newNumEvents >= N)
      break

    numEvents <- newNumEvents
  }

  return(index)
}

.makeBootstrapLenIndx <- function(s, indx.seq, N)
  sample(indx.seq, N, replace=TRUE)

sensitivitySGD <- function(z, s, d, y, beta0, beta1, phi, Pi, psi, tau,
                           time.points, 
                           selection, trigger, groupings,
                           ci, ci.method=c("analytic", "bootstrap"),
                           na.rm=FALSE, N.boot=100L, N.events=NULL,
                           oneSidedTest=FALSE, twoSidedTest=TRUE,
                           inCore=TRUE, verbose=getOption("verbose"),
                           colsPerFile=1000L) {
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

  ci.method <- sort(unique(match.arg(ci.method, several.ok=TRUE)))
  
  ErrMsg <- character(0L)
  if(isSlaveMode) {
    ## Running in boot strap mode
    
  } else {
    ## Not running a boot strap mode
    ## Run error checks on variables.
    ErrMsg <- NULL
    ErrMsg <- c(.CheckSelection(selection, s),
                .CheckGroupings(groupings),
                .CheckTrigger(trigger, d),
                .CheckLength(z=z, s=s, d=d, y=y),
                .CheckZ(z, groupings, na.rm),
                .CheckS(s, na.rm=na.rm),
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

    s <- s == selection
    d <- d == trigger
    
    z <- z == groupings[2L]

  }
  
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

  if(all(z0.s1 == FALSE) || all(z1.s1 == FALSE))
    if(isSlaveMode) {
      return(list(SCE = logical(0)))
    } else {
      stop("No events occured in one or more of the treatment arms")
    }

  if(!missing(psi) && !is.null(psi)) {
    Pi <-
      ifelse(abs(psi) < sqrt(.Machine$double.eps), p0*p1,
             -(sqrt((p1^2-2*p0*p1+p0^2)*exp(2*psi)+p1^2
                    +exp(psi)
                    *(-2*p1^2+2*p1-2*p0^2+2*p0)
                    +(2*p0-2)*p1+p0^2-2*p0+1)
               +p1+exp(psi)*(-p1-p0)+p0-1)
             /(2*exp(psi)-2))

    phi <- Pi/p1
  } else if(!missing(phi) && !is.null(phi)) {
    Pi <- p1*phi
    psi <- log((p1 * phi^2 + (1 - p0 - p1)*phi)/
               (p1 * phi^2 - (p1 + p0)* phi + p0))
  } else {
    psi <- log(Pi * (1 - p1 - p0 + Pi)/(p1 - Pi)/(p0 - Pi))
    phi <- Pi/p1
  }

  ## summary survfit of length of time til event for group 0 and group 1.
##   temp <- summary(survfit(Surv(y[z0.s1],d[z0.s1])~1L))
##   t0 <- temp$time
##   F0 <- 1L - temp$surv
   
##   temp <- summary(survfit(Surv(y[z1.s1],d[z1.s1])~1L))
##   t1 <- temp$time
##   F1 <- 1L - temp$surv

  ## Alt implementation
  ## survfit of length of time til event for group 0 and group 1.
  temp <- with(survfit(Surv(y[s], d[s]) ~ z[s], se.fit=FALSE), {
    who <- n.event > 0
    strata <- rep.int(seq_along(strata), times=strata)[who]
    time <- time[who]
    F <- 1L - surv[who]
    data.frame(strata=strata, time=time, F=F)
  })
  
  t0 <- temp$time[temp$strata==1]
  F0 <- temp$F[temp$strata==1]
  t1 <- temp$time[temp$strata==2]
  F1 <- temp$F[temp$strata==2]

  len.t0 <- length(t0)
  len.t1 <- length(t1)

  ## total length of phi vector
  len.total <- 4L + len.t0 + len.t1

  if(len.t0 == 0L || len.t1 == 0) {
    if(isSlaveMode) {
      return(list(SCE = logical(0)))
    } else {
      stop("No times occured in one or more of the treatment arms")
    }
  }

  dF0 <- diff(c(0L,F0,1L))
  dF1 <- diff(c(0L,F1,1L))

  coeffs <- .calc.coeff(Pi=Pi, p0=p0, p1=p1, beta0=beta0, beta1=beta1,
                       dF0=dF0, dF1=dF1, t0=t0, t1=t1, tau=tau,
                       time.points=time.points)

  SCE.dim <- c(length(time.points), length(beta0), length(beta1), length(psi))
  SCE.length <- prod(SCE.dim)
  SCE.dimnames <- list(time.points=as.character(time.points),
                       beta0=as.character(beta0), beta1=as.character(beta1),
                       psi=as.character(psi))
  SCE <- array(numeric(SCE.length), dim=SCE.dim, dimnames=SCE.dimnames)
  
  ## iterate across all betas and Pi values
  for(Pi.coeff in coeffs) {
    for(beta0.coeff in Pi.coeff$beta0.coeff) {
      for(beta1.coeff in Pi.coeff$beta1.coeff) {
        SCE[,beta0.coeff$i,beta1.coeff$i, Pi.coeff$i] <- beta0.coeff$SCE - beta1.coeff$SCE
      }
    }
  }
  
  
  if(isSlaveMode)
    return(list(SCE = SCE))

  cdfs <- list(alphahat0=sapply(coeffs, FUN=function(coeff) sapply(coeff$beta0.coeff, FUN=function(b) b$alphahat)),
               beta0=beta0,
               alphahat1=sapply(coeffs, FUN=function(coeff) sapply(coeff$beta1.coeff, FUN=function(b) b$alphahat)),
               beta1=beta1)

  z.seq <- seq_len(N)

  if(twoSidedTest) {
    if(ci < 0.5)
      ci.probs <- c(ci, 1L) - ci/2L
    else
      ci.probs <- c(0L, ci) + (1-ci)/2
  } else {
    ci.probs <- NULL
  }

  if(oneSidedTest) {
    ci.probs <- c(ci.probs, ci)
  }

  if("analytic" %in% ci.method) {
    stop("Analytic method is not currently implemented")
  }
  if("bootstrap" %in% ci.method) {
    current.fun <- sys.function()

    N.boot <- as.integer(N.boot)

    if(is.null(N.events)) {
      nVal <- N
      mkBsIndex <- .makeBootstrapLenIndx
    } else {
      nVal <- as.integer(N.events)
      mkBsIndex <- .makeBootstrapEvntIndx
    }

    bootCalc <- function(i, z.seq, nVal, beta0, beta1, psi, tau, time.points,
                         current.fun, verbose) {
      samp <- mkBsIndex(s, indx.seq=z.seq, N=nVal)
      ans <- as.vector(current.fun(z=z[samp], s=s[samp], d=d[samp],
                                   y=y[samp],
                                   beta0=beta0, beta1=beta1, psi=psi,
                                   tau=tau,
                                   selection=TRUE,
                                   groupings=c(FALSE,TRUE),
                                   trigger=TRUE,
                                   time.points=time.points,
                                   ci=NULL, ci.method=NULL)$SCE)
      if(verbose) cat(".")
      return(ans)
    }

    if(inCore) {
      vals <- apply(do.call(rbind, lapply(integer(N.boot), FUN=bootCalc,
                                          z.seq=z.seq, nVal=nVal,
                                          beta0=beta0, beta1=beta1,
                                          psi=psi, tau=tau,
                                          time.points=time.points,
                                          current.fun=current.fun,
                                          verbose=verbose)),
                    2L,
                    FUN=function(x) return(c(var(x), quantile(x, probs=ci.probs))))

      SCE.boot <- list(SCE.var = vals[1L,,drop=FALSE],
                       SCE.ci = vals[-1L,,drop=FALSE])
    } else {
      colsPerFile <- as.integer(colsPerFile)

      tmpfile <- tempfile()

      fieldWidth <- nchar(sprintf("%+a", .Machine$double.xmax))
      recordWidth <- fieldWidth + 1L
      lineWidth <- recordWidth*SCE.length
      fieldFmt <- sprintf("%%+%da%s", fieldWidth,
                          rep(c(" ", "\n"), times=c(SCE.length - 1L, 1L)))
      
      lapply(logical(N.boot),
             FUN=function(...) cat(sprintf(fieldFmt, bootCalc(...)), sep="", file=tmpfile,
               append=TRUE),
             z.seq=z.seq, nVal=nVal,
             beta0=beta0, beta1=beta1,
             psi=psi, tau=tau,
             time.points=time.points,
             current.fun=current.fun,
             verbose=verbose)
      if(verbose) cat('\n')

      
      needFiles <- (SCE.length %/% colsPerFile)
      remainder <- SCE.length %% colsPerFile

      filesNCols <- rep(c(colsPerFile, remainder), times=c(needFiles, remainder > 0L))
      readWidths <- filesNCols*recordWidth
      outFilenames <- sprintf("%s_split_%0*d", tmpfile,
                              nchar(length(readWidths)),
                              seq_along(readWidths))
      inConn <- file(tmpfile, open="r")
      on.exit(close(inConn))

      cullColumns <- function(offset, readWidth, outFilename, inConn,
                              lineWidth, N.boot) {
        outConn <- file(outFilename, open='w')
        on.exit(close(outConn))
        seek(inConn, where=offset, origin="start")

        lapply(logical(N.boot),
               FUN = function(j, inConn, outConn, readWidth, lineWidth) {
                 section <- readChar(inConn, readWidth)

                 if(substr(section, readWidth, readWidth) == " ")
                   substr(section, readWidth, readWidth) <- "\n"

                 writeChar(section, outConn, readWidth, eos=NULL)

                 seek(inConn, where=lineWidth - readWidth, origin="current")
                 return(NULL)
               }, inConn=inConn, outConn=outConn, readWidth=readWidth, lineWidth=lineWidth)

        cat("*")
        return(NULL)
      }
      
      mapply(FUN=cullColumns,
             offset=cumsum(c(0L, readWidths[-1L])),
             readWidth=readWidths,
             outFilename=outFilenames,
             MoreArgs=list(inConn=inConn, lineWidth=lineWidth, N.boot=N.boot))
      if(verbose) cat("\n")
      
      SCE.boot <- do.call(cbind, mapply(MoreArgs=list(ci.probs=ci.probs),
                                        FUN=function(filename, numCols, ci.probs) {
                                          dat <- scan(filename,
                                                      what=rep(list(double(0L)), numCols),
                                                      quiet=TRUE)
                                          ans <- sapply(dat,
                                                        FUN=function(column) c(var(column), quantile(column, ci.probs, names=TRUE)))
                                          
                                          if(verbose) cat('.')

                                          return(ans)
                                        },
                                        filename=outFilenames, numCols=filesNCols))
      if(verbose) cat('\n')
    }

    dim(SCE.boot$SCE.var) <- SCE.dim
    dimnames(SCE.boot$SCE.var) <- SCE.dimnames

    dim(SCE.boot$SCE.ci) <- c(nrow(SCE.boot$SCE.ci), SCE.dim)
    dimnames(SCE.boot$SCE.ci) <- c(list(ci.probs), SCE.dimnames)
    
  }

  ans <- structure(c(list(SCE=SCE), cdfs, SCE.boot))

  return(ans)
}
