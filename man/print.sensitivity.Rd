\name{print.sensitivity}
\alias{print.sensitivity}
\alias{print.sensitivity.0d}
\alias{print.sensitivity.1d}

\title{
  prints the results of calls to the sensitivity functions.
}
\description{
  Print the prints sensitivityPStrat objects in a visually understandable way.
}

\usage{
\method{print}{sensitivity.0d}(x, ...)
\method{print}{sensitivity.1d}(x, ...)
}
\arguments{
  \item{x}{
    sensitivity object
  }
  \item{\dots}{
    arguments passed to other print methods
  }
}

\author{
  Charles Dupont \cr
  Department of Biostatistics\cr
  Vanderbilt University\cr
}

\seealso{
  \code{\link{print.default}}  
}
\examples{
data(vaccine.trial)

print(with(vaccine.trial,
          sensitivityJR(z=treatment,s=hiv.outcome,y=logVL,
                    beta0=c(-1,-.75,-.5,-.25,0,.25,.5,.75,1),
                    beta1=c(-1,-.75,-.5,-.25,0,.25,.5,.75,1),
                    phi=c(0.95,0.90), selection="infected",
                    groupings=c("placebo","vaccine"),
                    N.boot=50)
         ))

}
\keyword{ print }
