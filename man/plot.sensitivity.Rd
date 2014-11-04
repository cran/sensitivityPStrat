\name{plot.sensitivity}
\alias{plot.sensitivity}
\alias{plot.sensitivity.1.0d}
\alias{plot.sensitivity.2.0d}
\alias{plot.sensitivity.1.1d}

\title{
  plots the results of calls to the sensitivity functions.
}
\description{
  Functions used to plot the objects created by the sensitivityPStrat
  family of functions.
}
\usage{
\method{plot}{sensitivity.1.0d}(x, xlim, ylim,
     xlab = expression(beta), ylab = "ACE",
     display = c("analytic", "bootstrap"),
     ci.select = 1,
     col = "black", line.col = col, point.col = col,
     analytic.col = "red", analytic.line.col = analytic.col,
     analytic.point.col = analytic.col,
     bootstrap.col = "green", bootstrap.line.col = bootstrap.col,
     bootstrap.point.col = bootstrap.col,
     panel.last = NULL, type = "l", ...)

\method{plot}{sensitivity.2.0d}(x, xlim, ylim, xlab = expression(beta[0]), ylab = expression(beta[1]),
     display = c("analytic", "bootstrap"), col = c(gray(.9), gray(1), gray(.8)),
     panel.last = NULL, ...)

\method{plot}{sensitivity.1.1d}(x, xlim, ylim,
     xlab = expression(beta), ylab = "SCE",
     t.point, display = c("analytic", "bootstrap"),
     col = "black", line.col = col, point.col = col,
     analytic.col = "red", analytic.line.col = analytic.col,
     analytic.point.col = analytic.col,
     bootstrap.col = "green", bootstrap.line.col = bootstrap.col,
     bootstrap.point.col = bootstrap.col,
     panel.last = NULL, type = "l", ...)
}
\arguments{
  \item{x}{
    sensitivity object
  }
  \item{t.point}{
    the time point at which data to create the plot.
  }
  \item{display}{
    character vector. Controls which confidence interval to use plot.
  }
  \item{ci.select}{
    integer vector or \samp{all}.  Selects the confidence interval to
     be ploted. If set to \samp{all} then all confidence intervals are
     plotted. Default value is 1.
   }
  \item{line.col}{
    the color all the lines should be.
  }
  \item{point.col}{
    the color all the infinity points should be.
  }
  \item{analytic.col}{
    vector; the color of all of the analytic confidence interval
    markings. Value are recycled if more confidence intervals are
    selected then given color values.
  }
  \item{analytic.line.col}{
    vector;the color of all of the analytic confidence interval lines.
    Value are recycled if more confidence intervals are selected then
    given color values.
  }
  \item{analytic.point.col}{
    vector; the color of all of the analytic confidence interval
    infinity points.  Value are recycled if more confidence intervals
    are selected then given color values.
  }
  \item{bootstrap.col}{
    vector; the color of all of the bootstrap confidence interval
     markings.  Value are recycled if more confidence intervals are
     selected then given color values.
  }
  \item{bootstrap.line.col}{
    vector; the color of all of the bootstrap confidence interval
     lines.  Value are recycled if more confidence intervals are
     selected then given color values.
  }
  \item{bootstrap.point.col}{
    vector; the color of all of the bootstrap confidence interval
     infinity points.  Value are recycled if more confidence intervals
     are selected then given color values.
  }
  \item{xlim, ylim, xlab, ylab, col, panel.last, type}{
    see \code{\link{plot.default}}
  }
  \item{\dots}{
    arguments passed to \code{\link{plot.default}}
  }
}
\author{
  Charles Dupont \cr
  Department of Biostatistics\cr
  Vanderbilt University\cr
}

\seealso{
  \code{\link{plot.default}}  
}
\examples{
data(vaccine.trial)

ansJR<-with(vaccine.trial,
          sensitivityJR(z=treatment,s=hiv.outcome,y=logVL,
                    beta0=c(-1,-.5,-.25,0,.25,.5,1),
                    beta1=c(-1,-.5,-.25,0,.25,.5,1),
                    phi=c(0.95,0.90), selection="infected",
                    groupings=c("placebo","vaccine"),
                    N.boot=50)
         )

plot(ansJR)

ans<-with(vaccine.trial,
          sensitivityGBH(z=treatment,s=hiv.outcome,y=logVL,
                    beta=c(-Inf,-1,0.75,-0.5,-0.25,0,.25,.5,.75,1,Inf),
                    selection="infected",
                    groupings=c("placebo","vaccine"),
                    empty.principal.stratum=c("not infected","infected"),
                    ci.method="bootstrap", ci=c(0.95, 0.9, 0.9),
                    ci.type=c("twoSided", "upper", "lower"),
                    custom.FUN=function(mu0, mu1, ...) mu1 - mu0,
                    N.boot=50, method=c("ACE", "T1", "T2"),
                    upperTest=TRUE, lowerTest=TRUE, twoSidedTest=TRUE)
         )

plot(ans, ci.select="all", bootstrap.col=c("red","green","blue"))
}
\keyword{ hplot }

