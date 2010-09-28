calc.coeff <- function(Pi, p0, p1, beta0, beta1, dF0, dF1, 
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

  coeffs <- mapply(FUN=calc.Pi.coeff,
                   Pi..p0=Pi..p0, Pi..p1=Pi..p1, i=seq.int(along.with=Pi),
                   MoreArgs=list(beta.tplus0=beta.tplus0, beta.tplus1=beta.tplus1,
                     q.seq.list=q.seq.list, r.seq.list=r.seq.list, dF0=dF0, dF1=dF1),
                   USE.NAMES = FALSE, SIMPLIFY=FALSE)
  return(coeffs)
}

calc.Pi.coeff <- function(Pi..p0, Pi..p1, i, beta.tplus0, beta.tplus1,
                          q.seq.list, r.seq.list, dF0, dF1, dV, dW, VmF0, WmF1,
                          p0, p1) {
  beta0.coeff <- lapply(beta.tplus0, q.list=q.seq.list, dF=dF0, Pi..p=Pi..p0,
                        FUN=calc.beta.coeff)
  beta1.coeff <- lapply(beta.tplus1, q.list=r.seq.list, dF=dF1, Pi..p=Pi..p1,
                        FUN=calc.beta.coeff)

  return(list(beta0.coeff=beta0.coeff, beta1.coeff=beta1.coeff, i=i))
}    

calc.beta.coeff <- function(beta.tplus, q.list, dV, dF, p, Pi..p, Pi.N) {

  calc.alphahat <- function(beta.tplus, dF, B) {
    ee <- function(a, beta.tplus, dF, B) {
      tmp <- sum((1L + exp(-a - beta.tplus))^(-1) * dF) - B
      tmp*tmp
    }

    range <- c(-100L,100L)
    while(TRUE) {
      alphahat <- optimize(ee, range, beta.tplus=beta.tplus, dF=dF,
                           B=B)$minimum

      limits <- range + c(-10L, 10L)
      if(alphahat < limits[2] && alphahat > limits[1])
        break
      
      warning("optimize overflow alphahat value invalid ", alphahat)

      if(alphahat > limits[2])
        range <- range + 75
      else if(alphahat < limits[1])
        range <- range - 75
      else
        stop("limit failure")      
    }

    return(alphahat)
  }

  alphahat <- calc.alphahat(beta.tplus$bt, dF=dF, B=Pi..p)

  if(alphahat > 90 || alphahat < -90) {
    warning("optimize overflow alphahat value invalid ", alphahat)
  }

  w <- exp(alphahat+beta.tplus$bt)/(1 + exp(alphahat+beta.tplus$bt))

  SCE <- sapply(q.list, w.dF=w*dF,
                FUN=function(q.seq, w.dF) sum(w.dF[q.seq])) / Pi..p

  return(list(i=beta.tplus$i, alphahat=alphahat, SCE=SCE))
}

makeBootstrapLenIndx <- function(s, indx.seq, N)
  sample(indx.seq, N, replace=TRUE)

sensitivitySGDFollowup <- function(z, s.star, v, d, y, beta0, beta1,
                                   phi, Pi, psi, tau, followup.time,
                                   time.points, ci=0.95,
                                   selection, groupings,
                                   empty.principal.stratum, trigger,
                                   ci.method=c("analytic", "bootstrap"),
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

  ci.method <- sort(unique(match.arg(ci.method, several.ok=TRUE)))
  
  if(na.rm == TRUE) {
    naIndex <- !(is.na(s.star) | is.na(z) | is.na(v) | (s.star & (is.na(d) | is.na(y))))

    z <- z[naIndex]
    v <- v[naIndex]
    s.star <- s.star[naIndex]
    d <- d[naIndex]
    y <- y[naIndex]
  }

  if(!is.null(ci)) {
    ErrMsg <- character(0)
    if(missing(empty.principal.stratum) || is.null(empty.principal.stratum) ||
       length(empty.principal.stratum) != 2L)
      ErrMsg <- c(ErrMsg, "'empty.principal.stratum' argument must be a two element vector")

    if(missing(selection) || is.null(selection) || length(selection) != 1L)
      ErrMsg <- c(ErrMsg, "'selection' argument must be a single element vector")
    else if(!(selection %in% empty.principal.stratum))
      ErrMsg <- c(ErrMsg, "'selection' value does not appear in specified levels of 's' from 'empty.principal.stratum'")

    if(missing(groupings) || is.null(groupings) || length(groupings) != 2L)
      ErrMsg <- c(ErrMsg, "'groupings' argument must be a two element vector")

    if(missing(trigger) || is.null(trigger) || length(trigger) != 1L)
      ErrMsg <- c(ErrMsg, "'trigger' argument must be a single element vector")

    if(any(is.na(z)))
      ErrMsg <- c(ErrMsg, "argument 'z' cannot contain any NA values")
    else if(length(unique(z)) > 2L)
      ErrMsg <- c(ErrMsg, "argument 'z' can contain at most 2 unique values")
    else if(sum(groupings %in% z) < 1L)
      ErrMsg <- c(ErrMsg, "None of the specified levels of 'z' from 'groupings' match supplied values of 'z'")

    
    if(any(is.na(s.star)))
      ErrMsg <- c(ErrMsg, "argument 's' cannot contain any NA values")
    else if(length(unique(s.star)) > 2L)
      ErrMsg <- c(ErrMsg, "argument 's' can contain at most 2 unique values")
    else if(sum(empty.principal.stratum %in% s.star) < 1L)
      ErrMsg <- c(ErrMsg, "None of the specified levels of 's' from 'empty.principal.stratum' match supplied values of 's'")

    if(any(s.star == selection & is.na(d)))
      ErrMsg <- c(ErrMsg, sprintf("argument 'd' cannont contain any NA value if the corisponding 's' is %s", selection))
    
    if(length(unique(d[!is.na(d)])) > 2L)
      ErrMsg <- c(ErrMsg, "argument 'd' can contain at most 2 unique non-NA values")
    else if(length(unique(d[!is.na(d)])) > 1L && all(unique(d[!is.na(d)]) != trigger))
      ErrMsg <- c(ErrMsg, "Value of 'trigger' does not match any value of 'd'")
    
    if(any(s.star == selection & is.na(y)))
      ErrMsg <- c(ErrMsg, sprintf("argument 'y' cannont contain any NA value if the corisponding 's' is %s", selection))
    
    if(length(ErrMsg) > 0L)
      stop(paste(ErrMsg, collapse="\n  "))
  }
  
  s.star <- s.star == selection
  d <- d == trigger

  s <- (v < followup.time) & s.star
  
  if(empty.principal.stratum[1L] == selection)
    z <- z == groupings[1L]
  else if(empty.principal.stratum[2L] == selection)
    z <- z == groupings[2L]

  ## N  - number subjects
  N <-length(z)

  z0.s1 <- !z & s
  z1.s1 <- z & s
  
  if(all(z0.s1 == FALSE) || all(z1.s1 == FALSE) ||
     all(d[z0.s1] == FALSE) || all(d[z1.s1] == FALSE))
    if(is.null(ci)) {
      return(list(SCE = logical(0)))
    } else {
      stop("No events occured in one or more of the treatment arms")
    }

  ## p0 - probiblity that subject in group 0 was selected
##  p0 <- 1-tail(x=summary(survfit(Surv(v[!z], s[!z])~1L, se.fit=FALSE))$surv,
##               n=1L)
  

  ## p1 - probablity that subject in group 1 was selected
##  p1 <- 1-tail(x=summary(survfit(Surv(v[z], s[z])~1L, se.fit=FALSE))$surv,
##               n=1L)

  ## alternate Implenation for p0 and p1
  temp <- with(survfit(Surv(v, s) ~ z, se.fit=FALSE),{
    who <- n.event > 0L
    data.frame(strata=rep.int(seq_along(strata), times=strata)[who],
               surv=surv[who])
  })

    
  p0 <- 1L - tail(x=temp$surv[temp$strata == 1L], n=1L)
  p1 <- 1L - tail(x=temp$surv[temp$strata == 2L], n=1L)
    
  ## Calc Pi because it will be needed later
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
##   temp <- summary(survfit(Surv(y[z0.s1],d[z0.s1])~1L, se.fit=FALSE))
##   t0 <- temp$time
##   F0 <- 1L - temp$surv
  
  
##   temp <- summary(survfit(Surv(y[z1.s1],d[z1.s1])~1L, se.fit=FALSE))
##   t1 <- temp$time
##   F1 <- 1L - temp$surv

  ## Alternate way to generate t0 F0 and t1 and F1
  temp <- with(survfit(Surv(y[s], d[s]) ~ z[s], se.fit=FALSE), {
    who <- n.event > 0L
    strata <- rep.int(seq_along(strata), times=strata)[who]
    time <- time[who]
    F <- 1L - surv[who]
    data.frame(strata=strata, time=time, F=F)
  })

  t0 <- temp$time[temp$strata == 1L]
  t1 <- temp$time[temp$strata == 2L]

  F0 <- temp$F[temp$strata == 1L]
  F1 <- temp$F[temp$strata == 2L]

  len.t0 <- length(t0)
  len.t1 <- length(t1)  

  if(len.t0 == 0L || len.t1 == 0) {
    if(isSlaveMode) {
      return(list(SCE = logical(0)))
    } else {
      stop("No times occured in one or more of the treatment arms")
    }
  }
  
  ## total length of phi vector
  len.total <- 4L + len.t0 + len.t1

  dF0 <- diff(c(0,F0,1))
  dF1 <- diff(c(0,F1,1))

  coeffs <- calc.coeff(Pi=Pi, p0=p0, p1=p1, beta0=beta0, beta1=beta1, dF0=dF0,
                       dF1=dF1, t0=t0, t1=t1, tau=tau, time.points=time.points)
  
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

  if(is.null(ci))
    return(list(SCE = SCE))

  ##  cdfs <- list(alphahat0=sapply(coeffs, FUN=function(coeff) sapply(coeff$beta0.coeff, FUN=function(b) b$alphahat)),
  ##               beta0=beta0, t0=t0, F0=F0,
  ##               alphahat1=sapply(coeffs, FUN=function(coeff) sapply(coeff$beta1.coeff, FUN=function(b) b$alphahat)),
  ##               beta1=beta1, t1=t1, F1=F1, N=N)
  cdfs <- NULL
  ansAttrs <- list()
  
  z.seq <- seq_len(N)

  if(twoSidedTest) {
    if(ci < 0.5)
      ci.probs <- c(ci, 1L) - ci/2L
    else
      ci.probs <- c(0L, ci) + (1L-ci)/2L
  } else {
    ci.probs <- NULL
  }

  if(oneSidedTest) {
    ci.probs <- c(ci.probs, ci)
  }

  if(ci.method == "analytic") {
    stop("Analytic method is not currently implemented")
  }
  if(ci.method == "bootstrap") {
    current.fun <- sys.function()

    N.boot <- as.integer(N.boot)

    if(is.null(N.events)) {
      nVal <- N
      mkBsIndex <- makeBootstrapLenIndx
    } else {
      nVal <- as.integer(N.events)
      mkBsIndex <- makeBootstrapEvntIndx
    }

    bootCalc <- function(i, z.seq, nVal, beta0, beta1, psi, tau, time.points,
                         current.fun, verbose) {
      samp <- mkBsIndex(s, indx.seq=z.seq, N=nVal)
      ans <- as.vector(current.fun(z=z[samp], s=s[samp], d=d[samp],
                                   y=y[samp], v=v[samp],
                                   beta0=beta0, beta1=beta1, psi=psi,
                                   tau=tau, followup.time=followup.time,
                                   selection=TRUE,
                                   groupings=c(FALSE,TRUE),
                                   empty.principal.stratum=c(FALSE,TRUE),
                                   trigger=TRUE,
                                   time.points=time.points,
                                   ci=NULL)$SCE)
      if(verbose) cat(".")
      return(ans)
    }

    if(inCore) {
      vals <- do.call(rbind, lapply(integer(N.boot), FUN=bootCalc,
                                    z.seq=z.seq, nVal=nVal,
                                    beta0=beta0, beta1=beta1,
                                    psi=psi, tau=tau,
                                    time.points=time.points,
                                    current.fun=current.fun,
                                    verbose=verbose))
      
      val.stats <- apply(vals, 2L,
                         FUN=function(x) return(c(var(x),
                           quantile(x, probs=ci.probs))))

      ansAttrs$bootReps <- N.boot
      ansAttrs$ActualBootReps <- nrow(vals)
      SCE.boot <- list(SCE.var = val.stats[1L,,drop=FALSE],
                       SCE.ci = val.stats[-1L,,drop=FALSE])
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
    dimnames(SCE.boot$SCE.ci) <- c(list(rownames(SCE.boot$SCE.ci)), SCE.dimnames)
    
  }

  
  ans <- c(list(SCE=SCE), cdfs, SCE.boot)
  attributes(ans) <- c(ansAttrs, list(names=names(ans)))

  return(ans)
}
