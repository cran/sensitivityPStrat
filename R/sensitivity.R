.calc.ecdf <- function(x) {
  n <- length(x)
  
  if (n < 1)
    stop("'x' must have 1 or more non-missing values")

  vals <- sort(unique(x))

  index <- match(x, vals)
  val.count = tabulate(index)
  return(list(F=cumsum(val.count)/n, vals=vals, index=index))
}

.foldUpperTri <- function(x) {
  xrow <- row(x)
  xcol <- col(x)
  xnrow <- nrow(x)
  upper <- xrow < xcol

  x[(xrow[upper] - 1)*xnrow + xcol[upper]] <- x[upper]
  return(x)
}

".sumCrossUpperTri<-" <- function(x, na.rm=TRUE, value) {
  xrow <- row(x)
  xcol <- col(x)
  indx <- xrow <= xcol
  
  x[indx] <- rowSums(value[xcol[indx],]*value[xrow[indx],], na.rm=na.rm)
  x
}

print.sensitivity1d <- function(x, ...) {
  labs <- attr(x, "parameters")
  if(all(c('s0','s1') %in% names(labs))) {
    cat("Empty Principle Stratum: ",
        paste("S(", labs$z0, ") = ",labs$s0, ", S(", labs$z1, ") = ",
              labs$s1, sep=''),
        "\n")
  }
  
  cat("ACE:\t", paste("E(Y(", labs$z1, ") - Y(", labs$z0,") | S(", labs$z0,
                      ") = S(", labs$z1, ") = ",labs$selected,")", sep=''),
      "\n")

  ACE <- x$ACE
  print(ACE)

  ci.method <- dimnames(x$ACE.ci)[['ci.method']]
  
  cat("\nACE confidence interval:\n")
  if("analytic" %in% ci.method) {
    cat("By analytic method\n")
    print(x$ACE.ci[,"analytic",])
  }

  if(all(c("bootstrap", "analytic") %in% ci.method))
     cat("\n")
  
  if("bootstrap" %in% ci.method) {
    cat("By bootstrap method, N = ", attr(x, 'N.boot'), "\n",sep='')
    print(x$ACE.ci[,"bootstrap",])
  }
}

print.sensitivity2d <- function(x, ...) {
  labs <- attr(x, "parameters")
  if(all(c('s0','s1') %in% names(labs))) {
    cat("Empty Principle Stratum: ",
        paste("S(", labs$z0, ") = ",labs$s0, ", S(", labs$z1, ") = ",
              labs$s1, sep=''),
        "\n")
  }
  
  cat("ACE:\t", paste("E(Y(", labs$z1, ") - Y(", labs$z0,") | S(", labs$z0,
                      ") = S(", labs$z1, ") = ",labs$selected,")", sep=''),
      "\n")

  ACE <- x$ACE
  print(ACE)

  ci.method <- dimnames(x$ACE.ci)[['ci.method']]

  cat("\nACE confidence interval:\n")
  if("analytic" %in% ci.method) {
    cat("By analytic method\n")
    print(x$ACE.ci[,"analytic",])
  }

  if(all(c("bootstrap", "analytic") %in% ci.method))
    cat("\n")
  
  if("bootstrap" %in% ci.method) {
    cat("By bootstrap method, N = ", attr(x, 'N.boot'), "\n",sep='')
    print(x$ACE.ci[,"bootstrap",])
  }

  invisible(NULL)
}


plot.sensitivity1d <- function(x, xlim, ylim, xlab=expression(beta), ylab='ACE',
                               display = c("analytic", "bootstrap"),
                               col='black', line.col=col, point.col=col,
                               analytic.col="red",
                               analytic.line.col=analytic.col,
                               analytic.point.col=analytic.col,
                               bootstrap.col="green",
                               bootstrap.line.col=bootstrap.col,
                               bootstrap.point.col=bootstrap.col,
                               panel.last=NULL,
                               type='l', ...) {
  display <- match.arg(display, several.ok=TRUE)

  sortIndx <- sort.list(x$beta)
  beta <- x$beta[sortIndx]
  ACE <- x$ACE[sortIndx]
  ACE.ci <- x$ACE.ci[sortIndx,,, drop=FALSE]

  finIndx <- is.finite(beta)
  infIndx <- is.infinite(beta)
  
  beta.fin <- beta[finIndx]
  ACE.fin <- ACE[finIndx]
  ACE.ci.fin <- ACE.ci[finIndx,,, drop=FALSE]

  beta.inf <- beta[infIndx]
  ACE.inf <- ACE[infIndx]
  ACE.ci.inf <- ACE.ci[infIndx,,, drop=FALSE]
  
  if(missing(ylim)) {
    ylim <- range(c(ACE, ACE.ci))
  }

  if(missing(xlim)) {
    xlim <- range(beta, finite=TRUE)
  }
  
  indx <- match(beta.inf, c(-Inf, Inf))

  plot.default(x=beta.fin, y=ACE.fin, xlim=xlim, ylim=ylim, type=type,
               ...)
  inf.x <- par('usr')[indx] + strwidth("m")/c(2,-2)[indx]
  points(x=inf.x, y=ACE.inf, pch=1, col=point.col)

  ##  doBootstrap
  if('analytic' %in% display && 'analytic' %in% colnames(x$ACE.var)) {
    for(i in seq_len(dim(ACE.ci.fin)[[2]])) {
      lines(x=beta.fin, y=ACE.ci.fin[,i,'analytic'], lty=2,
            col=analytic.line.col)
    }
    points(x=rep.int(inf.x, times=dim(ACE.ci.fin)[[2]]),
           y=ACE.ci.inf[,indx, "analytic"], pch=3, col=analytic.point.col)
  }

  if('bootstrap' %in% colnames(x$ACE.var) && 'bootstrap' %in% display) {
    for(i in seq_len(dim(ACE.ci.fin)[[2]])) {
      lines(x=beta.fin, y=ACE.ci.fin[,i,"bootstrap"], lty=2,
            col=bootstrap.line.col)
    }
    points(x=rep(inf.x, times=dim(ACE.ci.fin)[[2]]),
           y=ACE.ci.inf[,indx,"bootstrap"], pch=3, col=bootstrap.point.col)
  }
}

plot.sensitivity2d <- function(x, xlim, ylim, xlab=expression(beta),
                                 t.point,
                                 ylab='SCE',
                                 display = c("analytic", "bootstrap"),
                                 col='black', line.col=col, point.col=col,
                                 analytic.col="red",
                                 analytic.line.col=analytic.col,
                                 analytic.point.col=analytic.col,
                                 bootstrap.col="green",
                                 bootstrap.line.col=bootstrap.col,
                                 bootstrap.point.col=bootstrap.col,
                                 panel.last=NULL,
                                 type='l', ...) {
  
  display <- match.arg(display, several.ok=TRUE)

  sortIndx <- sort.list(x$beta)
  beta <- x$beta[sortIndx]
  SCE <- x$SCE[sortIndx]
  SCE.ci <- x$SCE.ci[t.point, sortIndx,,, drop=FALSE]

  finIndx <- is.finite(beta)
  infIndx <- is.infinite(beta)
  
  beta.fin <- beta[finIndx]
  SCE.fin <- SCE[finIndx]
  SCE.ci.fin <- SCE.ci[t.point, finIndx,,, drop=FALSE]

  beta.inf <- beta[infIndx]
  SCE.inf <- SCE[infIndx]
  SCE.ci.inf <- SCE.ci[t.point, infIndx,,, drop=FALSE]
  
  if(missing(ylim)) {
    ylim <- range(c(SCE, SCE.ci))
  }

  if(missing(xlim)) {
    xlim <- range(beta, finite=TRUE)
  }
  
  indx <- match(beta.inf, c(-Inf, Inf))

  plot.default(x=beta.fin, y=SCE.fin, xlim=xlim, ylim=ylim, type=type,
               ...)
  inf.x <- par('usr')[indx] + strwidth("m")/c(2,-2)[indx]
  points(x=inf.x, y=SCE.inf, pch=1, col=point.col)

  ##  doBootstrap
  if('analytic' %in% display && 'analytic' %in% colnames(x$SCE.var)) {
    lines(x=beta.fin, y=SCE.ci.fin[t.point,,1,'analytic'], lty=2, col=analytic.line.col)
    lines(x=beta.fin, y=SCE.ci.fin[t.point,,2,'analytic'], lty=2, col=analytic.line.col)
    points(x=rep.int(inf.x, times=2), y=SCE.ci.inf[t.point,indx,, "analytic"], pch=3, col=analytic.point.col)
  }

  if('bootstrap' %in% colnames(x$SCE.var) && 'bootstrap' %in% display) {
    lines(x=beta.fin, y=SCE.ci.fin[t.point,,1,"bootstrap"], lty=2, col=bootstrap.line.col)
    lines(x=beta.fin, y=SCE.ci.fin[t.point,,2,"bootstrap"], lty=2, col=bootstrap.line.col)
    points(x=rep(inf.x, times=2), y=SCE.ci.inf[t.point,indx,,"bootstrap"], pch=3, col=bootstrap.point.col)
  }
}
