sensitivityJR <- function(z, s, y, beta0, beta1, phi, Pi, psi,
                          selection, groupings, ci=0.95,
                          ci.method=c("analytic", "bootstrap"), na.rm=FALSE,
                          N.boot=100, oneSidedTest = FALSE,
                          twoSidedTest = TRUE, verbose=getOption("verbose"))
{
  w.calc <- function(alpha, beta, y)
    1/(1 + exp(-alpha - beta*y))
  
  alpha.est <- function(alpha, beta, y, dF, C)
    (sum(w.calc(alpha=alpha, beta=beta, y=y) * dF) - C)^2


  calc.coefs <- function(y, beta, dF, RR) {
    coefs <- vector(length(beta), "list")

    for(i in seq_along(beta)) {
      alphahat <- alpha.est(y=y, beta=beta[i], dF=dF, RR=RR)
      w <- w.calc(alpha=alphahat, beta=beta[i], y=y)
      coefs[[i]] <- list(alphahat=alphahat, w=w)
    }

    return(coefs)
  }

  if(!missing(ci.method) && is.null(ci.method))
    isSlaveMode <- TRUE
  else
    isSlaveMode <- FALSE

  ci.method <- sort(unique(match.arg(ci.method, several.ok=TRUE)))

  if(!isSlaveMode) {

    if(na.rm == TRUE) {
      naIndex <- !(is.na(s) | is.na(z))

      z <- z[naIndex]
      s <- s[naIndex]
      y <- y[naIndex]
    }
    
    ErrMsg <- character(0)
    if(missing(selection) || is.null(selection))
      ErrMsg <- c(ErrMsg, "'selection' argument must be a single element vector")
    else if(!selection %in% s)
      ErrMsg <- c(ErrMsg, "value of 'selection' is not included in 's'")
    
    if(missing(groupings) || is.null(groupings))
      ErrMsg <- c(ErrMsg, "'groupings' argument must be a two element vector")
    else if(!all(groupings %in% z))
      ErrMsg <- c(ErrMsg, "Specified levels of 'z' from 'groupings' do not match supplied values of 'z'")

    if(any(is.na(z)))
      ErrMsg <- c(ErrMsg, "argument 'z' cannot contain any NA values")
    else if(length(unique(z)) != 2)
      ErrMsg <- c(ErrMsg, "argument 'z' can only contain 2 unique values")

    
    if(any(is.na(s)))
      ErrMsg <- c(ErrMsg, "argument 's' cannot contain any NA values")
    else if(length(unique(s)) > 2)
      ErrMsg <- c(ErrMsg, "argument 's' can contain at most 2 unique values")

    ParamConflict <- c(!missing(psi) && !is.null(psi),
                       !missing(Pi) && !is.null(Pi),
                       !missing(phi) && !is.null(phi))

    if(sum(ParamConflict) > 1)
      ErrMsg <- c(ErrMsg,
                  "only one of 'psi', 'Pi', or 'phi' can be specified")
    else if(sum(ParamConflict) < 1)
      ErrMsg <- c(ErrMsg,
                  "one of 'psi', 'Pi', or 'phi' must be specified")

    if(length(ErrMsg) > 0)
      stop(paste(ErrMsg, collapse="\n  "))

    s <- s == selection

    z <- z == groupings[2]
    
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
  }
  
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

  RR <- p1/p0
  VE <- 1 - RR

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

  ## Calc psi for later use
    
  Fn0 <- ecdf(y0)
  y0.uniq <- knots(Fn0)
  F0 <- Fn0(y0.uniq)
  dF0 <- diff(c(0, F0))
  
  Fn1 <- ecdf(y1)
  y1.uniq <- knots(Fn1)
  F1 <- Fn1(y1.uniq)
  dF1 <- diff(c(0, F1))

#  print(Pi)
  q0c <- matrix(quantile(y0, probs=c(Pi/p0, 1-Pi/p0)), ncol=2)
  q1c <- matrix(quantile(y1, probs=c(Pi/p1, 1-Pi/p1)), ncol=2)
  
  ACE.dim <- c(length(beta0), length(beta1), length(Pi))
  ACE.length <- prod(ACE.dim)
  ACE.dimnames <- list(beta0=format(beta0, trim=TRUE),
                       beta1=format(beta1, trim=TRUE),
                       Pi = format(Pi, trim=TRUE))

  ACE <- array(numeric(ACE.length), dim=ACE.dim, dimnames=ACE.dimnames)

  FnAs0.dim <- ACE.dim[-2L]
  FnAs0.length <- prod(FnAs0.dim)
  FnAs0.dimnames <- ACE.dimnames[-2L]

  mu0 <- alphahat0 <- array(numeric(FnAs0.length), dim=FnAs0.dim,
                            dimnames=FnAs0.dimnames)

  FnAs1.dim <- ACE.dim[-1L]
  FnAs0.length <- prod(FnAs0.dim)
  FnAs1.dimnames <- ACE.dimnames[-1L]

  mu1 <- alphahat1 <- array(numeric(FnAs0.length), dim=FnAs1.dim,
                            dimnames=FnAs1.dimnames)

  if(!isSlaveMode) {
    FnAs0 <- funArray(vector(mode='list', length=prod(FnAs0.dim)),
                      dim=FnAs0.dim,
                      dimnames=FnAs0.dimnames)

    FnAs1 <- funArray(vector(mode='list', length=prod(FnAs1.dim)),
                      dim=FnAs1.dim,
                      dimnames=FnAs1.dimnames)
  }

  for(i in seq_along(Pi)) {
    if(Pi[i] == 0) {
      ACE[,,i] <- NA
      alphahat0[,i] <- NA
      alphahat1[,i] <- NA

      next
    }
    
    a0 <- Pi[i]/p0
    a1 <- Pi[i]/p1
    
    for(j in seq_along(beta0)) {
      if(is.finite(beta0[j])) {
        alphahat0[j,i] <- optimize(alpha.est, c(-100,100), beta=beta0[j],
                                   y=y0.uniq, dF=dF0, C=a0)$minimum
        w0 <- w.calc(alpha=alphahat0[j,i], beta=beta0[j], y=y0.uniq)

        dFas0 <- dF0*w0/a0
        if(!isSlaveMode)
          Fas0 <- cumsum(dFas0)
        
        mu0[j,i] <- sum(y0.uniq * dFas0)
      } else {
        if(beta0[j] == Inf) {
          Fas0 <- ifelse(y0.uniq <= q0c[i,1L] & F0 < a0, F0/a0, 1)
        } else if(beta0[j] == -Inf) {
          Fas0 <- ifelse(y0.uniq >= q0c[i,2L], (F0 - (1 - a0))/a0, 0)
        } else {
          stop("Invalid beta0 value ", beta0[j])
        }

        alphahat0[j,i] <- NA
        mu0[j,i] <- sum(y0.uniq * diff(c(0, Fas0)))
      }

      if(!isSlaveMode)
        FnAs0[j,i] <- stepfun(y0.uniq, c(0, Fas0))

    }
    
    for(j in seq_along(beta1)) {
      if(is.finite(beta1[j])) {
        alphahat1[j,i] <- optimize(alpha.est, c(-100,100), beta=beta1[j],
                                   y=y1.uniq, dF=dF1, C=a1)$minimum
        
        w1 <- w.calc(alpha=alphahat1[j,i], beta=beta1[j], y=y1.uniq)

        dFas1 <- dF1*w1/a1

        if(!isSlaveMode)
          Fas1 <- cumsum(dFas1)
        
        mu1[j,i] <- sum(y1.uniq * dFas1)
      } else if(is.infinite(beta1[j])) {
        if(beta1[j] == Inf) {
          Fas1 <- ifelse(y1.uniq <= q1c[i,1L] & F1 < a1, F1/a1, 1L)
        } else if(beta1[j] == -Inf) {
          Fas1 <- ifelse(y1.uniq >= q1c[i,2L], (F1 - a1)/(a1), 0L)
        }

        alphahat1[j,i] <- NA
        mu1[j,i] <- sum(y1.uniq * diff(c(0L, Fas1)))
      } else {
        stop("Invalid beta1 value ", beta1[j])
      }

      if(!isSlaveMode)
        FnAs1[j,i] <- stepfun(y1.uniq, c(0, Fas1))
    }

    ACE[,,i] <- outer(mu0[,i], mu1[,i], function(mu0, mu1) mu1 - mu0)
  }
  
  if(isSlaveMode)
    return(list(ACE=ACE))

  cdfs<-list(beta0=beta0, alphahat0=alphahat0, Fas0=FnAs0,
             beta1=beta1, alphahat1=alphahat1, Fas1=FnAs1,
             Pi=Pi, phi=phi, psi=psi)

  if(is.null(ci))
    return(c(list(ACE=ACE), cdfs))
  
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

  ci.probs <- sort(unique(ci.probs))
  ci.probsLen <- length(ci.probs)
  
  ACE.ci.dim <- c(ACE.dim, length(ci.method), ci.probsLen)
  ACE.ci.length <- prod(ACE.ci.dim)
  ACE.ci.dimnames <- c(ACE.dimnames,
                       list(ci.method=ci.method,
                            ci.probs=sprintf("%s%%",
                                             as.character(ci.probs*100))))
  
  ACE.ci <- array(numeric(ACE.ci.length), dim=ACE.ci.dim,
                  dimnames=ACE.ci.dimnames)

  ACE.var.dim <- c(ACE.dim, length(ci.method))
  ACE.var.length <- prod(ACE.var.dim)
  ACE.var.dimnames <- c(ACE.dimnames, list(ci.method=ci.method))

  ACE.var <- array(numeric(ACE.var.length), dim=ACE.var.dim,
                   dimnames=ACE.var.dimnames)

  ## run bootstrap method
  if(any(ci.method == "analytic")) {
    if(verbose)
      cat("running Analytic")
    
    Omega <- matrix(nrow=6,ncol=6)
    for(k in seq_along(Pi)) {
      for(i in seq_along(beta0)) {
        for(j in seq_along(beta1)) {
          U <- rbind((1-z)*(p0 - s),
                     (z)*(p1 - s),
                     (1-z)*s*(1/(1+exp(-alphahat0[i,k] - beta0[i]*y)) - Pi[k]/p0),
                     (z)*s*(1/(1+exp(-alphahat1[j,k] - beta1[j]*y)) - Pi[k]/p1),
                     (1-z)*s*(mu0[i,k] - y*p0/Pi[k]/(1+exp(-alphahat0[i,k] - beta0[i]*y))),
                     (z)*s*(mu1[j,k] - y*p1/Pi[k]/(1+exp(-alphahat1[j,k] - beta1[j]*y))))
          sumCrossUpperTri(Omega) <- U
          Omega <- Omega / N
          Omega[lower.tri(Omega)] <- t(Omega)[lower.tri(Omega)]

          Gamma <- matrix(colSums(cbind(1-z,
                                        0,
                                        s*(1-z)*Pi[k]/p0^2,
                                        0,
                                        -s*y*(1-z)/Pi[k]/(exp(-beta0[i]*y-alphahat0[i,k])+1),
                                        0,

                                        
                                        0,
                                        z,
                                        0,
                                        s*z*Pi[k]/p1^2,
                                        0,
                                        -s*y*z/Pi[k]/(exp(-beta1[j]*y-alphahat1[j,k])+1),

                                        
                                        0,
                                        0,
                                        s*exp(-beta0[i]*y-alphahat0[i,k])*(1-z)/(exp(-beta0[i]*y-alphahat0[i,k])+1)^2,
                                        0,
                                        -p0*s*y*exp(-beta0[i]*y-alphahat0[i,k])*(1-z)/Pi[k]/(exp(-beta0[i]*y-alphahat0[i,k])+1)^2,
                                        0,

                                        
                                        0,
                                        0,
                                        0,
                                        s*exp(-beta1[j]*y-alphahat1[j,k])*z/(exp(-beta1[j]*y-alphahat1[j,k])+1)^2,
                                        0,
                                        -p1*s*y*exp(-beta1[j]*y-alphahat1[j,k])*z/Pi[k]/(exp(-beta1[j]*y-alphahat1[j,k])+1)^2,

                                        
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
                                        0,
                                        s*z), na.rm=TRUE),
                          nrow=6,byrow=TRUE) / N

          IGamma <- solve(Gamma)
          ## vartheta <- tcrossprod(IGamma %*% Omega, IGamma) / N
          vartheta <- (t(IGamma) %*% Omega %*% IGamma) / N
          ## cat(all.equal(vartheta, alt), '\n',sep='')
          ACE.var[i, j, k, "analytic"] <- vartheta[5, 5] + vartheta[6, 6] - 2 * vartheta[5, 6]

          if(verbose) cat(".")
        }
      }
    }

    calculateCi <- function(i, norm, ACE, sqrt.ACE.var) {
      ACE[i] + norm * sqrt.ACE.var[i]
    }
    
    ACE.ci[,,,"analytic",] <- outer(seq_along(ACE), qnorm(ci.probs),
                                   FUN=calculateCi, ACE=ACE,
                                   sqrt.ACE.var=sqrt(ACE.var[,,,'analytic']))
    
    if(verbose) cat("\n")
  }

  if(any(ci.method == "bootstrap")) {
    if(verbose) cat("running Bootstrap")
    ACE.list <- array(dim=c(ACE.dim,N.boot))
    
    for(i in seq_len(N.boot)) {
      new.index <- sample(seq_len(N), N, replace=TRUE)
      ACE.list[,,,i] <- Recall(z=z[new.index],s=s[new.index],y=y[new.index],
                               beta0=beta0, beta1=beta1, psi=psi, ci=NULL,
                               ci.method=NULL)$ACE
      if(verbose) cat(".")
    }

    for(k in seq_along(Pi)) {
      for(i in seq_along(beta0)) {
        for(j in seq_along(beta1)) {
          ACE.ci[i,j,k,'bootstrap',] <- quantile(ACE.list[i,j,k,], probs=ci.probs)
          ACE.var[i,j,k,'bootstrap'] <- var(ACE.list[i,j,k,])
        }
      }
    }
    if(verbose) cat("\n")
  }

  return(structure(c(list(ACE=ACE, ACE.ci=ACE.ci, ACE.var=ACE.var),
                     cdfs), class="sensitivity3d", N.boot=N.boot,
                   parameters=list(z0=groupings[1], z1=groupings[2],
                     selected=selection)))
}

