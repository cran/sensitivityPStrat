\name{sensitivitySGL}
\alias{sensitivitySGL}
\title{principal stratification sensitivity analysis with time to event data}
\description{
  principal stratification sensitivity analysis with time to event data
  using method described in SGL.
}
\usage{
sensitivitySGL(z, s, d, y, beta, tau, time.points, ci=0.95, selection,
               groupings, empty.principle.stratum, trigger,
               ci.method = c("analytic", "bootstrap"), na.rm = FALSE,
               N.boot = 100L, oneSidedTest = FALSE, twoSidedTest = TRUE,
               verbose = getOption("verbose"))
}
\arguments{
  \item{z}{
    vector containing the grouping values of the the records.
  }
  \item{s}{
    vector describing whether an record is selected.
  }
  \item{d}{
    vector describing whether an event has occured for this record. Can
    be \code{NA} for unselected records
  }
  \item{y}{
    vector describing the length after the event occured until observed
    effect.  Can be \code{NA} for unselected records or records where
    there was no event.
  }
  \item{beta}{
    Vector of beta values of intrest.
  }
  \item{tau}{
    max length of time between event and observed effect.
  }
  \item{time.points}{
    vector of time points where a SCE value is desired.
  }
  \item{ci}{
    numeric.  confindence interval.
  }
  \item{selection}{
    The value of \code{s} in where a record is selected.
  }
  \item{groupings}{
    Vector of two values, first being the value of \code{z} which
    delineates the control group, the last being the value of \code{z}
    which delineates the treatment group.
  }
  \item{empty.principle.stratum}{
    vector of two. Describing the \code{s} values which select the empty
    principle stratum.  If empty.principle.stratum=c(s0,s1), then
    stratum defined by S(groupings[1])==s0 and S(groupings[2])==s1 is
    the empty stratum.  In this example s1 and s2 refer to the two
    possible values of s.
  }
  \item{trigger}{
    logical.  The value of \code{d} which denotes an event occurance.
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

  \item{SCE}{
    array; Calculated values of SCE for all combinations of the values from
    \code{phi}/\code{Pi}/\code{psi}, and
    \code{time.points}.  Array dimentions are
    \code{length(time.points)}, \code{length(psi)}.
  }
  \item{beta}{
    beta values used.
  }
  \item{psi}{
    phi values used.
  }
  \item{SCE.ci}{
    array; Confindence interval of the SCE value.  Confidance interval
    determined by \code{quantile} if using  \code{ci.method}
    \dQuote{bootstrap}.  Other wise calculated using analytic
    varance. Array dimentions the same as element \code{SCE}.
  }
  \item{SCE.var}{
    array; Variance of the SCE value.  Array dimentions the same as
    element \code{SCE}.
  }
}
\references{
  Shepherd BE, Gilbert PB, Lumley T (2007), "Sensitivity analyses
  comparing time-to-event outcomes existing only in a subset selected
  postrandomization," Journal of the American Statistical Association
  102, 573-582.
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
  \code{\link{sensitivitySGD}}, \code{\link{sensitivityJR}}
  \code{\link{Surv}}
}

\examples{
#sensitivitySGL(z=bladder$rx, s=bladder$enum < 5,
#               d=bladder$event, y=bladder$stop,
#               beta=c(-1,0,1), 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{design }
