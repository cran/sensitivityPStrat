\name{sensitivitySGL}
\alias{sensitivitySGL}
\title{principal stratification sensitivity analysis with time to event data}
\description{
  Principal stratification sensitivity analysis with time to event data
  using the method described by Shepherd, Gilbert, and Lumley (2007).
}
\usage{
sensitivitySGL(z, s, d, y, v, beta, tau, time.points, selection, trigger,
               groupings, empty.principal.stratum, followup.time,
               ci=0.95, ci.method = c("analytic", "bootstrap"),
               ci.type="twoSided", custom.FUN = NULL, na.rm = FALSE,
               N.boot = 100L, interval = c(-100, 100),
               upperTest = FALSE, lowerTest = FALSE, twoSidedTest = TRUE,
               verbose = getOption("verbose"), isSlaveMode = FALSE)
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
  \item{beta}{
    vector; values of the sensitivity parameter \eqn{\beta}.  \code{Inf} and
    \code{-Inf} are acceptable.
  }
  \item{tau}{
    maximum observed follow-up time after selection.  Selection weights
    are constant for \eqn{t>\code{tau}}.
  }
  \item{time.points}{
    vector; time points, \var{t}, at which \eqn{SCE(\var{t})} will be estimated.
  }
  \item{selection}{
    The value of \code{s} indicating selection.
  }
  \item{trigger}{
    logical; the value of \code{d} that denotes the post-selection event.
  }
  \item{groupings}{
    Vector of two elements \code{c(\var{g0},\var{g1})}, the first element \var{g0} being the
    value of \code{z} the delineates the first group, the last element
    \var{g1} being the value of \code{z} which delineates the second group.
  }
  \item{empty.principal.stratum}{
    vector of two elements \code{c(\var{s0},\var{s1})}; describes the \code{s}
    values that select the empty principal stratum.  If
    \code{empty.principal.stratum=c(\var{s0},\var{s1})}, then stratum defined by
    \eqn{S(\var{g0})=\var{s0}} and \eqn{S(\var{g1})=\var{s1}} is the empty stratum.  In this example \var{s0}
    and \var{s1} refer to the two possible values of \code{s}. (Note: method only
    works if \eqn{\var{s0} \ne \var{s1}}{\var{s0} != \var{s1}}).
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
    \dQuote{bootstrap}. 
  }
  \item{ci.type}{
    character vector; type of confidence interval that the corresponding
    \code{ci} element is referring to.  Can be \dQuote{upper},
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
  \item{interval}{
    numeric vector of length 2. Controls the range limits used to
    by optimize to estimate \eqn{\alpha}.
  }
  \item{lowerTest}{
    logical;  Return the lower one sided p-value for SCE. Defaults
    to \code{FALSE}
  }
  \item{upperTest}{
    logical;  Return the upper one sided p-value for SCE. Defaults
    to \code{FALSE}
  }
  \item{twoSidedTest}{
    logical;  Return a two sided p-value for SCE. Defaults
    to \code{TRUE}
  }
  \item{verbose}{
    logical;  prints dots when bootstrapping to show that something is
    happening.  Bootstrapping can take a long time.
  }
  \item{isSlaveMode}{
    logical.  Internal Use only. Used in recursion.
  }
}
\details{
  Performs a sensitivity analysis estimating the \dQuote{survival causal
    effect} among those who would have been selected regardless of
  treatment assignment (SCE).  The method assumes no interference (i.e.,
  potential outcomes of all subjects are unaffected by treatment
  assignment of other subjects), ignorable (i.e., random) treatment
  assignment, monotonicity (i.e., one of the principal strata is
  empty), and independent censoring (i.e., time from selection to
  event is independent of time from selection until censoring).  SCE is
  then identified by assuming a value of the sensitivity parameter \eqn{\beta},
  where \eqn{e^{\beta}}{exp(\beta)} has an odds ratio interpretation (see help for
  \code{\link{sensitivityGBH}}).  Given selection in one treatment arm, the
  probability of selection if in the other treatment arm is assumed to be constant for
  for \eqn{T(\code{z})>\code{tau}}.  
  
  SCE is computed at user specified time points.
  
  Specifying \code{beta}=\code{-Inf} or
  \code{beta}=\code{Inf} estimates the bounds for SCE.
}
\value{
  object of class \code{sensitivity2d}

  \item{SCE}{
    \eqn{SCE(\var{t})=Pr(T(\var{g0})\le\var{t}|S(\var{g0})=S(\var{g1})=\code{selection}) - Pr(T(\var{g1})\le\var{t}|S(\var{g0})=S(\var{g1})=\code{selection})}{SCE(\var{t})=Pr(T(\var{g0})<=\var{t}|S(\var{g0})=S(\var{g1})=\code{selection}) - Pr(T(\var{g1})<=\var{t}|S(\var{g0})=S(\var{g1})=\code{selection})}.
    Array of the estimated SCE at all \code{time.points} for specified
    beta values.  Array dimensions are  \code{length(time.points)} by \code{length(beta)}.
  }
  \item{SCE.ci}{
    array; confidence interval of SCE
    determined by \code{quantile} if using \code{ci.method} includes 
    \dQuote{bootstrap}.  Otherwise calculated using analytic
    variance with large sample normal approximation. Array dimensions the same as element \code{SCE}.
  }
  \item{SCE.var}{
    array; estimated variance of SCE.  Array dimensions the same as
    element \code{SCE}.
  }
  \item{ci.map}{
    list; mapping of confidence interval to quantile probability.  Use
    numbers contained within as indices to the \code{SCE.ci} element.
  }
  \item{beta}{
    vector of user-specified \eqn{\beta} values
  }
  \item{alphahat}{
    vector of estimated values of \eqn{\alpha}
  }
  \item{y0}{
    vector of unique event times in the first group.
  }
  \item{Fas0}{
    matrix of estimated empirical distribution function values for
    \code{y0} in the first group in the always selected principal
    stratum. \eqn{Pr(Y(\var{g0})\le\var{y0}|S(\var{g0})=S(\var{g1})=\code{selection}; \beta)}{Pr(Y(\var{g0}) <= \var{y0}|S(\var{g0})=S(\var{g1})=\code{selection}; \beta)}
  }
  \item{y1}{
    vector of unique event times in the second group.
  }
  \item{Fas1}{
    matrix of estimated empirical distribution function values for
    \code{y1} in the second group in the always selected principal
    stratum.  \eqn{Pr(Y(\var{g1})\le\var{y1}|S(\var{g0})=S(\var{g1})=\code{selection}; \beta)}{Pr(Y(\var{g1}) <= \var{y1}|S(\var{g0})=S(\var{g1})=\code{selection}; \beta)}
  }
}
\references{
  Shepherd BE, Gilbert PB, Lumley T (2007), \dQuote{Sensitivity analyses
  comparing time-to-event outcomes existing only in a subset selected
  postrandomization,} Journal of the American Statistical Association
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
  \code{\link{sensitivityGBH}}, \code{\link{sensitivityHHS}}, \code{\link{sensitivitySGD}}, 
  \code{\link{Surv}}
}

\examples{

data(vaccine.trial)
sens.time<-with(vaccine.trial,
                sensitivitySGL(z=treatment, s=hiv.outcome, y=followup.yearsART,
                          d=ARTinitiation, beta=c(.25, 0,-.25), tau=3,
                          time.points=c(2,3), selection="infected",
                          trigger="initiated ART", groupings=c("placebo","vaccine"),
                          empty.principal.stratum=c("not infected","infected"),
                          N.boot=50, interval=c(-200,200),
                          upperTest=TRUE, lowerTest=TRUE, twoSidedTest=TRUE)
               )
sens.time

sens.time2<-with(vaccine.trial,
                sensitivitySGL(z=treatment, s=hiv.outcome, y=followup.yearsART,
                          d=ARTinitiation, beta=c(.25, 0,-.25), tau=3,
                          time.points=c(2,3), selection="infected",
                          trigger="initiated ART", groupings=c("placebo","vaccine"),
                          empty.principal.stratum=c("not infected","infected"),
                          custom.FUN=function(Fas0,Fas1,time.points,
                ...) { Fas0(time.points) - Fas1(time.points) },
                          N.boot=50, interval=c(-200,200))
               )
sens.time2

sens.time3<-with(vaccine.trial,
                sensitivitySGL(z=treatment, s=hiv.outcome, y=followup.yearsART,
                          d=ARTinitiation, beta=c(-Inf,.25,0,-.25,Inf),
                          tau=3, time.points=c(2,3), selection="infected",
                          trigger="initiated ART", groupings=c("placebo","vaccine"),
                          empty.principal.stratum=c("not infected","infected"),
                          custom.FUN=function(Fas0,Fas1,time.points,
                ...) { Fas0(time.points) - Fas1(time.points) },
                          N.boot=50, interval=c(-200,200))
               )
sens.time3

}
