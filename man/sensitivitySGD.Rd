\name{sensitivitySGD}
\alias{sensitivitySGD}
\title{
  principal stratification sensitivity analysis with time to event data
  relaxing monotonicity assumption.
}
\description{
  Principal stratification sensitivity analysis with time to event data
  relaxing monotonicity as described by Shepherd, Gilbert, and Dupont
  (in press).
}
\usage{
sensitivitySGD(z, s, d, y, v, beta0, beta1, phi, Pi, psi, tau,
               time.points, selection, trigger, groupings,
               followup.time,
               ci=0.95, ci.method = c("bootstrap", "analytic"),
               ci.type="twoSided", custom.FUN = NULL, na.rm = FALSE,
               N.boot = 100L, N.events = NULL, interval = c(-100, 100),
               upperTest = FALSE, lowerTest = FALSE, twoSidedTest=TRUE,
               inCore = TRUE,verbose = getOption("verbose"),
               colsPerFile = 1000L, isSlaveMode = FALSE)
}
\arguments{
  \item{z}{
    vector; contains the grouping values (e.g., treatment assignment)
    for each record.
  }
  \item{s}{
    vector; indicates whether a record is selected.
  }
  \item{d}{
    vector; indicates whether a post-selection event has occurred. Can
    be \code{NA} for unselected records.
  }
  \item{y}{
    vector; the length of time from selection until event (\code{d}) or
    censoring.
    Can be \code{NA} for unselected records.
  }
  \item{v}{
    numeric vector; the length of time from randomization until
    selection or censoring.
  }
  \item{beta0}{
    numeric vector; values of the sensitivity parameter \eqn{\beta} linking outcome in group
    \var{g0} with \code{selection} if assigned group \var{g1}.
  }
  \item{beta1}{
    numeric vector; values of the sensitivity parameter \eqn{\beta} linking outcome in group
    \var{g1} with \code{selection} if assigned group \var{g0}.
  }
  \item{phi, Pi, psi}{
    vectors; sensitivity parameters specifying the joint distribution of \eqn{S(\var{g0})},
    \eqn{S(\var{g1})}.  Only one of the three parameters should be
    specified. \code{psi} is the log-odds ratio of selection.
    \code{Pi} is the probability of being in the always selected
    principal stratum (\eqn{Pr(S(\var{g0}) = S(\var{g1}) = selected)}).
    \code{phi} is the probability of selection in group \var{g0} given
    selection in group \var{g1} (\eqn{Pr(S(\var{g0}) = 1|S(\var{g1}) = 1)}).
  }
  \item{tau}{
    maximum observed follow-up time after selection.  Selection weights
    are constant for \eqn{\var{t}>\code{tau}}.
  }
  \item{time.points}{
    vector; time points, \var{t}, at which \eqn{SCE(\var{t})} will be estimated.
  }
  \item{selection}{
    The value of \code{s} indicating selection.
  }
  \item{trigger}{
    The value of \code{d} that denotes the post-selection event.
  }
  \item{groupings}{
    Vector of two elements \code{c(\var{g0},\var{g1})}, the first element \var{g0} being the
    value of \code{z} the delineates the first group, the last element
    \var{g1} being the value of \code{z} which delineates the second group.
  }
  \item{followup.time}{
    numeric value; cut-off point for \code{v} after which records are
    lost to censoring.
  }
  \item{ci}{
    numeric vector; confidence interval level, defaults to \code{0.95}.
  }
  \item{ci.method}{
    character;  method by which the confidence interval and
    variance are calculated.  Can be \dQuote{analytic} or
    \dQuote{bootstrap}. Currently only works for \dQuote{bootstrap}.
  }
  \item{ci.type}{
    character vector; type of confidence interval that the corresponding
    \code{ci} element is refering to.  Can be \dQuote{upper},
    \dQuote{lower}, or \dQuote{twoSided}.  Defaults to
    \code{"twoSided"}.
  }    
  \item{custom.FUN}{
    function; function to calculate custom result. \code{Fas0}, \code{Fas1},
    \code{time.points}, \code{p0}, \code{p1} are available to be used as
    arguments in the custom function.  The custom function must return a
    vector of elements that is the same length as \code{time.points}.
  }
  \item{na.rm}{
    logical; indicates whether records that are invalid due to \code{NA}
    values should be removed from the data set.
  }
  \item{N.boot}{
    integer; number of bootstrap repetitions that will be run when 
    \code{ci.method} includes \dQuote{bootstrap}.
  }
  \item{N.events}{
    integer; number of selection-events (S) for each bootstrap
    replication when doing selection-event based bootstrapping.
  }
  \item{interval}{
    numeric vector of length 2. Controls the range limits used to
    by \code{\link{optimize}} to estimate \eqn{\alpha}.
  }
  \item{lowerTest}{
    logical.  Return the lower one sided p-value for SCE. Defaults
    to \code{FALSE}
  }
  \item{upperTest}{
    logical.  Return the upper one sided p-value for SCE. Defaults
    to \code{FALSE}
  }
  \item{twoSidedTest}{
    logical.  Return a two sided p-value for SCE. Defaults
    to \code{TRUE}
  }
  \item{verbose}{
    logical;  prints dots when bootstrapping to show that something is
    happening.  Bootstrapping can take a long time.
  }
  \item{inCore}{
    logical; running in memory if \code{TRUE}, running with scratch files if
    \code{FALSE}.  Default is \code{TRUE}.  For large data analysis, the
    user may want to switch this to \code{FALSE} to allow for processing
    on data sets larger than can fit in memory.
  }
  \item{colsPerFile}{
    integer; number of columns of the scratch file to process in each
    pass (e.g., 100 columns).
  }
  \item{isSlaveMode}{
    logical.  Internal Use only. Used in recursion.
  }
}
\details{
  Performs a sensitivity analysis estimating the \dQuote{survival causal
    effect} among those who would have been selected regardless of
  treatment assignment (SCE) without assuming monotonicity (i.e., that
    one of the principal stratum is empty).  The method assumes no interference (i.e.,
  potential outcomes of all subjects are unaffected by treatment
  assignment of other subjects), ignorable (i.e., random) treatment
  assignment, and independent censoring (i.e., time from selection to
  event is independent of time from selection until censoring).  SCE is
  identified by assuming values for the sensitivity parameters \code{beta0},
    \code{beta1}, and one of the parameters \code{phi}, \code{psi}, or \code{Pi}.  The sensitivity
    parameters \code{beta0} and \code{beta1} have a log-odds ratio interpretation (see
    help for \code{\link{sensitivityGBH}}).  Given selection in one treatment arm, the
  probability of selection if in the other treatment arm is assumed to be constant for
  for \eqn{T(\code{z})>\code{tau}}.
  
  Only one of the parameters \code{phi}, \code{psi}, or \code{Pi} should be specified as all
  depend on each other.  \code{psi} is unrestrained taking any value on the real
  line.  The other parameters, \code{phi} and \code{Pi} have constraints and there
  will be estimation problems if these parameters are set at values
  outside the of their range of acceptable values based on the observed
  data.  See Shepherd, Gilbert, Dupont (in press) for more details.
}
\value{
  object of class \code{sensitivity3d}

  \item{SCE}{
    array; Calculated values of SCE for all combinations of the values from
    \code{beta0}, \code{beta1}, \code{phi}/\code{Pi}/\code{psi}, and
    \code{time.points}.  Array dimensions are
    \code{length(time.points)}, \code{length(beta0)},
    \code{length(beta1)}, \code{length(psi)}.
  }
  \item{SCE.ci}{
    array; Confidence interval of the SCE value.  Confidence interval
    determined by \code{\link{quantile}} if using \code{ci.method}
    \dQuote{bootstrap}.  Otherwise calculated using analytic
    variance with large sample normal approximation. Array dimensions the same as element \code{SCE}.
  }
  \item{SCE.var}{
    array; estimated variance of SCE.  Array dimensions the same as
    element \code{SCE}.
  }
  \item{beta0}{
    vector; \eqn{\beta} values used for first group.
  }
  \item{beta1}{
    vector; \eqn{\beta} values used for second group.
  }
  \item{psi}{
    vector; \eqn{\psi} values used.
  }
  \item{Pi}{
    vector; \eqn{Pi} values used.
  }
  \item{psi}{
    vector; \eqn{psi} values used.
  }
  \item{ci.map}{
    list; mapping of confidence interval to quantile probability.  Use
    numbers contained within as indices to the \code{SCE.ci} element.
  }
}
\references{
  Shepherd BE, Gilbert PB, and Dupont CT, \dQuote{Sensitivity analyses
  comparing time-to-event outcomes only existing in a subset selected
  postrandomization and relaxing monotonicity,} Biometrics, in press.
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
  \code{\link{sensitivitySGL}}, \code{\link{sensitivityJR}},
  \code{\link[survival]{Surv}}
}

\examples{

data(vaccine.trial)
sens.analysis<-with(vaccine.trial,
                sensitivitySGD(z=treatment, s=hiv.outcome, y=followup.yearsART,
                          d=ARTinitiation, beta0=c(0,-.25,-.5),
                          beta1=c(0, -.25, -.5), phi=c(0.95, 0.90), tau=3,
                          time.points=c(2,3), selection="infected",
                          trigger="initiated ART",
                          groupings=c("placebo","vaccine"), ci=.95,
                          ci.method="bootstrap", N.boot=100)
               )
sens.analysis


sens.analysis2<-with(vaccine.trial,
                sensitivitySGD(z=treatment, s=hiv.outcome, y=followup.yearsART,
                          d=ARTinitiation, beta0=c(0,-.25,-.5),
                          beta1=c(0, -.25, -.5), phi=c(0.95, 0.90), tau=3,
                          time.points=c(2,3), selection="infected",
                          trigger="initiated ART",
                          groupings=c("placebo","vaccine"), ci=.95,
                          custom.FUN=function(Fas0,Fas1,...,time.points) {
                            Fas0(time.points) - Fas1(time.points)
                          },
                          ci.method="bootstrap", N.boot=100)
               )
sens.analysis2

}

