\name{sensitivityJR}
\alias{sensitivityJR}
\title{
  Sensitivity analysis using the JR method.
}
\description{
  Performs a sensitivity analysis using the method described in JR.
}
\usage{
sensitivityJR(z, s, y, beta0, beta1, phi, Pi, psi, selection,
              groupings, ci = 0.95, ci.method = c("analytic",
                                                  "bootstrap"),
              na.rm = FALSE, N.boot = 100, oneSidedTest = FALSE,
              twoSidedTest = TRUE, verbose=getOption("verbose"))
}
\arguments{
  \item{z}{
    vector containing the grouping values of the the records.
  }
  \item{s}{
    vector describing whether an record is selected.
  }
  \item{y}{
    vector describing the length after the event occured until observed
    effect.  Can be \code{NA} for unselected records or records where
    there was no event.
  }
  \item{beta0, beta1}{
    vector beta values for group 0 and group 1.  \code{Inf} and
    \code{-Inf} are acceptable values.
  }
  \item{phi}{
    Don't Know
  }
  \item{Pi}{
    Don't Know
  }
  \item{psi}{
    Don't Know
  }
  \item{selection}{
    The value of \code{s} in where a record is selected.
  }
  \item{groupings}{
    Vector of two values, first being the value of \code{z} which
    delineates the control group, the last being the value of \code{z}
    which delineates the treatment group.
  }
  \item{ci}{
    numeric vector;  confindence interval.
  }
  \item{ci.method}{
    character.  Method by which the convidence interval and
    variance are calculated.  Can be either \dQuote{analytic} or
    \dQuote{bootstrap}. 
  }
  \item{na.rm}{
    logical. Should \code{NA} be removed from the data set.
  }
  \item{N.boot}{
    integer. Number of bootstrap repitions that will be run when running with
    \code{ci.method} equal to \dQuote{bootstrap}.
  }
  \item{oneSidedTest}{
    logical.  Return a one sided confidence interval.
  }
  \item{twoSidedTest}{
    logical.  Return a two sided confidence interval.
  }
  \item{verbose}{
    logical. Print dots when boot straping to show that something is happening.
  }
}
\value{
  object of class \code{sensitivity3d}

  \item{ACE}{
    array of the calculated ACE values for the data set and
    bound selection.  Array has dimentions of \code{length(beta0)},
    \code{length(beta1)}, \code{length(phi)}.
  }
  \item{beta0}{
    beta0 values used.
  }
  \item{beta1}{
    beta1 values used.
  }
  \item{phi}{
    phi values used.
  }
  \item{ACE.ci}{
    array; confidance interval determined by \code{quantile} if using
    \code{ci.method} \dQuote{bootstrap}.  Other wise calculated using
    analytic varance.  Has same dimentions as \code{ACE} element.
  }
  \item{ACE.var}{
    array; Variance of the analysis.   Has same dimentions as \code{ACE}
    element.
  }
}
\references{
  Jemiai Y (2005), "Semiparametric Methods for Inferring Treatment
  Effects on Outcomes Defined Only if a Post-Randomization Event
  Occurs," unpublished doctoral dissertation under the supervision of
  A. Rotnitzky, Harvard School of Public Health, Dept. of Biostatistics.
}
\author{
  Bryan E. Shepherd \cr
  Department of Biostatistics\cr
  Vanderbilt University\cr

  
  Charles Dupont \cr
  Department of Biostatistics\cr
  Vanderbilt University\cr
}
\seealso{
  \code{\link[stats]{quantile}}, \code{\link[stats]{var}},
  \code{\link[sensitivityPStrat]{calc.v}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{design}
