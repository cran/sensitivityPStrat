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
    character vector. Controlls which confidence interval to use plot.
  }
  \item{line.col}{
    the color all the lines should be.
  }
  \item{point.col}{
    the color all the infinity points should be.
  }
  \item{analytic.col}{
    the color of all of the analytic confidence interval markings
  }
  \item{analytic.line.col}{
    the color of all of the analytic confidence interval lines    
  }
  \item{analytic.point.col}{
    the color of all of the analytic confidence interval infinity points
  }
  \item{bootstrap.col}{
    the color of all of the bootstrap confidence interval markings
  }
  \item{bootstrap.line.col}{
    the color of all of the bootstrap confidence interval lines    
  }
  \item{bootstrap.point.col}{
    the color of all of the bootstrap confidence interval infinity points
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
                    beta0=c(-1,-.75,-.5,-.25,0,.25,.5,.75,1),
                    beta1=c(-1,-.75,-.5,-.25,0,.25,.5,.75,1),
                    phi=c(0.95,0.90,0.80), selection="infected",
                    groupings=c("placebo","vaccine"),
                    N.boot=100)
         )

plot(ansJR)

}
\keyword{ hplot }

