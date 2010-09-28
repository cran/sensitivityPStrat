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
sensitivitySGD(z, s, d, y, beta0, beta1, phi, Pi, psi, tau,
               time.points, selection, trigger, groupings,
               ci, ci.method = c("analytic", "bootstrap"),
               na.rm = FALSE, N.boot = 100L, N.events = NULL,
               oneSidedTest = FALSE, twoSidedTest = TRUE, inCore = TRUE,
               verbose = getOption("verbose"), colsPerFile = 1000L)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{beta0}{
    vector; values of the sensitivity parameter linking outcome in group
    g0 with selection if assigned group g1.
  }
  \item{beta1}{
    vector; values of the sensitivity parameter linking outcome in group
    g1 with selection if assigned group g0.
  }
  \item{phi, Pi, psi}{
    vectors; sensitivity parameters specifying the joint distribution of S(g0),
    S(g1).  Only one of the three parameters should be
    specified. \code{psi} is the log-odds ratio of selection.
    \code{Pi} is the probability of being in the always selected
    principal stratum (Pr(S(g0)=S(g1)=selected)).
    \code{phi} is the probability of selection in group g0 given
    selection in group g1 (Pr(S(g0)=1|S(g1)=1)).
  }
  \item{tau}{
    maximum observed follow-up time after selection.  Selection weights
    are constant for t>\code{tau}.
  }
  \item{time.points}{
    vector; time points, t, at which SCE(t) will be estimated.
  }
  \item{selection}{
    The value of \code{s} indicating selection.
  }
  \item{trigger}{
    logical; the value of \code{d} that denotes the post-selection event.
  }
  \item{groupings}{
    Vector of two elements \code{c(g0,g1)}, the first element \code{g0} being the
    value of \code{z} the delineates the first group, the last element
    \code{g1} being the value of \code{z} which delineates the second group.
  }
  \item{ci}{
    numeric vector; confidence interval level, defaults to 0.95.
  }
  \item{ci.method}{
    character;  method by which the confidence interval and
    variance are calculated.  Can be \dQuote{analytic} or
    \dQuote{bootstrap}. Currently only works for \dQuote{bootstrap}.
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
  \item{oneSidedTest}{
    logical.  Return a one sided confidence interval for
    SCE(t). Defaults to \code{FALSE}.
  }
  \item{twoSidedTest}{
    logical.  Return a two sided confidence interval for SCE(t).
    Defaults to \code{TRUE}.
  }
  \item{verbose}{
    logical;  prints dots when bootstrapping to show that something is happening.
  }
  \item{inCore}{
    logical; running in memory if TRUE, running with scratch files if
    FALSE.  Default is TRUE.  For large data analysis, the user may want
    to switch this to FALSE to make computations faster.
  }
  \item{colsPerFile}{
    integer; number of columns of the scratch file to process in each
    pass (e.g., 100 columns).
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
  identified by assuming values for the sensitivity parameters beta0,
    beta1, and one of the parameters phi, psi, or Pi.  The sensitivity
    parameters beta0 and beta1 have a log-odds ratio interpretation (see
    help for sensitivityGBH).  Given selection in one treatment arm, the
  probability of selection if in the other treatment arm is assumed to be constant for
  for T(z)>tau.  
  
  Only one of the parameters phi, psi, or Pi should be specified as all
  depend on each other.  psi is unrestrained taking any value on the real
  line.  The other parameters, psi and Pi have constraints and there
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
  \item{beta0}{
    beta0 values used.
  }
  \item{beta1}{
    beta1 values used.
  }
  \item{psi}{
    psi values used.
  }
  \item{SCE.ci}{
    array; Confidence interval of the SCE value.  Confidence interval
    determined by \code{quantile} if using  \code{ci.method}
    \dQuote{bootstrap}.  Otherwise calculated using analytic
    variance with large sample normal approximation. Array dimensions the same as element \code{SCE}.
  }
  \item{SCE.var}{
    array; estimated variance of SCE.  Array dimensions the same as
    element \code{SCE}.
  }
}
\references{
  Shepherd BE, Gilbert PB, and Dupont CT, "Sensitivity analyses
  comparing time-to-event outcomes only existing in a subset selected
  postrandomization and relaxing monotonicity," Biometrics, in press.
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
                          ci.method="bootstrap", N.boot=500)
               )
sens.analysis
}

