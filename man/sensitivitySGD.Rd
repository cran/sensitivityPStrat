\name{sensitivitySGD}
\alias{sensitivitySGD}
\title{
  principal stratification sensitivity analysis using the HHS method.
}
\description{
  Performs a sensitivity analysis on time to event data using method
  described in SGL.
}
\usage{
sensitivitySGD(z, s, d, y, beta0, beta1, phi, Pi, psi, tau, time.points,
               ci, selection, groupings,
               trigger, ci.method = c("analytic", "bootstrap"),
               na.rm = FALSE, N.boot = 100L, N.events = NULL,
               oneSidedTest = FALSE, twoSidedTest = TRUE, inCore = TRUE,
               verbose = getOption("verbose"), colsPerFile = 1000L)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{beta0}{
    Vector of beta0 values of intrest.
  }
  \item{beta1}{
    Vector of beta1 values of intrest.
  }
  \item{phi, Pi, psi}{
    Don't Know
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
  \item{N.events}{
    integer. Number of events to select when doing event based bootstrapping.
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
  \item{inCore}{
    Should this analysis be run using scratch files? (much much slower
    how ever very large analysis are allowed
  }
  \item{colsPerFile}{
    How many columns of the scratch file should be process in each pass.
  }
}
\value{
  object of class \code{sensitivity3d}

  \item{SCE}{
    array; Calculated values of SCE for all combinations of the values from
    \code{beta0}, \code{beta1}, \code{phi}/\code{Pi}/\code{psi}, and
    \code{time.points}.  Array dimentions are
    \code{length(time.points)}, \code{length(beta0)},
    \code{length(beta1)}, \code{length(psi)}.
  }
  \item{beta0}{
    beta0 values used.
  }
  \item{beta1}{
    beta1 values used.
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
  Shepherd BE, Gilbert PB, and Dupont CT, "Sensitivity analyses
  comparing time-to-event outcomes only existing in a subset selected
  postrandomization and relaxing monotonicity," (submitted paper).
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
  \code{\link[survival]{Surv}}, \code{\link[survival]{survfit}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{design}
\keyword{survival}% __ONLY ONE__ keyword per line
