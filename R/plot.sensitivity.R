plot.sensitivity.1.0d <- function(x, xlim, ylim,
                                  xlab=expression(beta), ylab='ACE',
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
               xlab=xlab, ylab=ylab, ...)
  
  inf.x <- par('usr')[indx] + strwidth("m")/c(2,-2)[indx]
  points(x=inf.x, y=ACE.inf, pch=1, col=point.col)

  ##  doBootstrap
  if('analytic' %in% display && 'analytic' %in% colnames(x$ACE.var)) {
    for(i in seq_len(dim(ACE.ci.fin)[[2]])) {
      lines(x=beta.fin, y=ACE.ci.fin[, i, 'analytic'], lty=2,
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

plot.sensitivity.2.0d <- function(x, xlim, ylim,
                                  xlab=expression(beta[0]),
                                  ylab=expression(beta[1]),
                                  display = c("analytic", "bootstrap"),
                                  col=c(gray(.9),gray(1),gray(.8)),
                                  panel.last=NULL, ...) {

  display <- match.arg(display, several.ok=FALSE)

  sortIndx0 <- sort.list(x$beta0)
  sortIndx1 <- sort.list(x$beta1)
  
  beta0 <- x$beta0[sortIndx0]
  beta1 <- x$beta1[sortIndx1]
  
  ACE <- x$ACE[sortIndx0,sortIndx1,]
  
  ACE.ci <- x$ACE.ci[sortIndx0,sortIndx1,,,, drop=FALSE]

  reject <- ifelse(ACE.ci[,,,2,display] < 0,
                  -1,
                  ifelse(ACE.ci[,,,1,display] > 0, 
                         1,
                         0))

  finIndx0 <- is.finite(beta0)
  infIndx0 <- is.infinite(beta0)

  finIndx1 <- is.finite(beta1)
  infIndx1 <- is.infinite(beta1)

  beta0.fin <- beta0[finIndx0]
  beta1.fin <- beta1[finIndx1]

  reject.fin <- reject[finIndx0,finIndx1,, drop=FALSE]

  beta0.inf <- beta0[infIndx0]
  beta1.inf <- beta1[infIndx1]

  reject.inf <- reject[infIndx0,infIndx1,,drop=FALSE]

  if(missing(ylim)) {
    ylim <- range(beta1, finite=TRUE)
  }

  if(missing(xlim)) {
    xlim <- range(beta0, finite=TRUE)
  }
  

  for(i in seq_len(dim(reject)[3])) {
    image(x=beta0.fin, y=beta1.fin, z=reject.fin[,,i], xlim=xlim, ylim=ylim,
          xlab=xlab, ylab=ylab, breaks=c(-1.5,-0.5,0.5,1.5), axes=FALSE,
          col=col,
          sub=bquote(.(as.symbol(names(dimnames(ACE))[3])) == .(format(as.numeric(dimnames(ACE)[[c(3,i)]]), digits=3))),
          ...)

    contour(beta0.fin,beta1.fin,reject.fin[,,i],add=TRUE, lty=3,
            levels=seq(from=-1, to=1, length.out=21),,labcex=0.8, axes=FALSE)
    axis(side=1,at=beta0,label=format(beta0, digits=1,trim=TRUE,drop0trailing=TRUE), line=NA)
    axis(side=2,at=beta1,label=format(beta1, digits=1,trim=TRUE,drop0trailing=TRUE), line=NA)
  }
  
#  inf.x <- par('usr')[indx] + strwidth("m")/c(2,-2)[indx]
#  points(x=inf.x, y=ACE.inf, pch=1, col=point.col)
}

plot.sensitivity.1.1d <- function(x, xlim, ylim, xlab=expression(beta), ylab="SCE",
                                 t.point,
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
  SCE.ci <- x$SCE.ci[sortIndx, t.point,,, drop=FALSE]

  finIndx <- is.finite(beta)
  infIndx <- is.infinite(beta)
  
  beta.fin <- beta[finIndx]
  SCE.fin <- SCE[finIndx]
  SCE.ci.fin <- SCE.ci[finIndx, t.point,,, drop=FALSE]

  beta.inf <- beta[infIndx]
  SCE.inf <- SCE[infIndx]
  SCE.ci.inf <- SCE.ci[infIndx, t.point,,, drop=FALSE]
  
  if(missing(ylim)) {
    ylim <- range(c(SCE, SCE.ci))
  }

  if(missing(xlim)) {
    xlim <- range(beta, finite=TRUE)
  }
  
  indx <- match(beta.inf, c(-Inf, Inf))

  plot.default(x=beta.fin, y=SCE.fin, xlim=xlim, ylim=ylim, type=type,
               xlab=xlab, ylab=ylab, ...)
  inf.x <- par('usr')[indx] + strwidth("m")/c(1,-1)[indx]
  points(x=inf.x, y=SCE.inf, pch=1, col=point.col)

  ##  doBootstrap
  if('analytic' %in% display && 'analytic' %in% colnames(x$SCE.var)) {
    for(i in seq_len(dim(SCE.ci.fin)[[3]])) {
      lines(x=beta.fin, y=SCE.ci.fin[, t.point, i,'analytic'], lty=2,
            col=analytic.line.col)
      points(x=inf.x, y=SCE.ci.inf[,t.point, i, "analytic"], pch=3, col=analytic.point.col)
    }
  }

  if('bootstrap' %in% colnames(x$SCE.var) && 'bootstrap' %in% display) {
    for(i in seq_len(dim(SCE.ci.fin)[[3]])) {
      lines(x=beta.fin, y=SCE.ci.fin[,t.point,i,'bootstrap'], lty=2,
            col=analytic.line.col)
    }
    points(x=rep(inf.x, times=2), y=SCE.ci.inf[indx, t.point,,"bootstrap"], pch=3, col=bootstrap.point.col)
  }
}
