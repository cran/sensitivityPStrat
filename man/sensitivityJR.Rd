\name{sensitivityJR}
\alias{sensitivityJR}
\title{
  Principal stratification sensitivity analysis relaxing the
  monotonicity assumption.
}
\description{
  Principal stratification sensitivity analysis relaxing monotonicity as
  described by Jemiai and Rotnitzky (2005) and implemented by Shepherd,
  Redman, and Ankerst (2008).
}
\usage{
sensitivityJR(z, s, y, beta0, beta1, phi, Pi, psi,
              selection, groupings,
              ci = 0.95, ci.method = c("analytic","bootstrap"),
              ci.type = "twoSided", custom.FUN=NULL, na.rm = FALSE,
              N.boot = 100, interval = c(-100, 100),
              upperTest = FALSE, lowerTest = FALSE, twoSidedTest = TRUE,
              verbose=getOption("verbose"), isSlaveMode = FALSE)
}
\arguments{
  \item{z}{
    vector; contains the grouping values (e.g., treatment assignment)
    for each record.
  }
  \item{s}{
    vector; indicates whether a record is selected.
  }
  \item{y}{
    vector; outcome values.  Can be \code{NA} for unselected records.
  }
  \item{beta0}{
    vector; values of the sensitivity parameter \var{\eqn{\beta}0} linking outcome in group
    \var{g0} with selection if assigned group \var{g1}.
  }
  \item{beta1}{
    vector; values of the sensitivity parameter \var{\eqn{\beta}1} linking outcome in group
    \var{g1} with selection if assigned group \var{g0}.
  }
  \item{phi, Pi, psi}{
    vector; sensitivity parameters specifying the joint distribution of
    \eqn{S(\var{g0})}, \eqn{S(\var{g1})}.  Only one of the three
    parameters should be specified. \code{psi} is the log-odds ratio of
    selection. \code{Pi} is the probability of being in the always
    selected principal stratum (\eqn{Pr(S(\var{g0}) = S(\var{g1}) =
      selected)}). \code{phi} is the probability of selection in group
    \var{g0} given selection in group \var{g1}
    (\eqn{Pr(S(\var{g0}) = 1|S(\var{g1}) = 1)}).
  }
  \item{selection}{
    The value of \code{s} indicating selection.
  }
  \item{groupings}{
    vector of two elements \code{c(\var{g0},\var{g1})}; describes to
    possible group values. The first element \var{g0} being the value of
    \code{z} the delineates the first group, the last element \var{g1}
    being the value of \code{z} which delineates the second group.
  }
  \item{ci}{
    numeric vector;  confidence interval value. Defaults to \code{0.95}
  }
  \item{ci.method}{
    character;  method by which the confidence interval and
    variance are calculated.  Can be \dQuote{analytic} or
    \dQuote{bootstrap}. Defaults to \code{c("analytic","bootstrap")}
  }
  \item{ci.type}{
    character vector; type of confidence interval that the corresponding
    \code{ci} element is referring to.  Can be \dQuote{upper},
    \dQuote{lower}, or \dQuote{twoSided}.  Defaults to
    \code{"twoSided"}.
  }
  \item{custom.FUN}{
    function; function to calculate custom result. \code{mu0},
    \code{mu1}, \code{p0}, \code{p1} are available to be used as
    arguments in the custom function, where
    \eqn{\code{mu0} = E(Y(\var{g0})|S(\var{g0}) = S(\var{g1}) = selected)},
    \eqn{\code{mu1} = E(Y(\var{g1})|S(\var{g0}) = S(\var{g1}) = selected)},
    \eqn{\code{p0} = P(S(\var{g0}) = selected)}, and
    \eqn{\code{p1} = P(S(\var{g1}) = selected)}.  The custom function
    must return a single value.
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
    numeric vector of length 2. Controls the range limits used
    by optimize to estimate \var{\eqn{\alpha}0} and \var{\eqn{\alpha}1}.
  }
  \item{lowerTest}{
    logical.  Return the lower one sided p-value for the ACE. Defaults
    to \code{FALSE}
  }
  \item{upperTest}{
    logical.  Return the upper one sided p-value for the ACE. Defaults
    to \code{FALSE}
  }
  \item{twoSidedTest}{
    logical.  Return a two sided p-value for the ACE. Defaults
    to \code{TRUE}
  }
  \item{verbose}{
    logical;  prints dots when bootstrapping to show that something is
    happening. Bootstrapping can take a long time.
  }
  \item{isSlaveMode}{
    logical.  Internal Use only. Used in recursion.
  }
}
\details{
  Performs a sensitivity analysis estimating the average causal effect
  among those who would have been selected regardless of treatment
  assignment (ACE) without assuming monotonicity (i.e., that one of the
  principal strata is empty).  The method assumes no interference (i.e.,
  potential outcomes of all subjects are unaffected by treatment
  assignment of other subjects) and ignorable (i.e., random) treatment
  assignment.  ACE is identified by assuming values for the sensitivity
  parameters \code{beta0}, \code{beta1}, and one of the parameters \code{phi}, \code{psi}, or \code{Pi}.
  The sensitivity parameters \code{beta0} and \code{beta1} have a log-odds ratio
  interpretation (see help for \code{\link{sensitivityGBH}}).
  
  Only one of the parameters \code{phi}, \code{psi}, or \code{Pi} should
  be specified as all depend on each other.  \code{psi} is unrestrained
  taking any value on the real line.  The other parameters, \code{psi}
  and \code{Pi} have constraints and there will be estimation problems
  if these parameters are set at values outside the of their range of
  acceptable values based on the observed data.  See Shepherd, Gilbert,
  Dupont (in press) for more details.
}
\value{
  object of class \code{sensitivity3d}

  \item{ACE}{
    array; estimated values of ACE for all combinations of \code{beta0},
    \code{beta1}, and \code{phi}, \code{Pi}, \code{psi}.  Array
    dimensions are \code{length(beta0)}, \code{length(beta1)},
    \code{length(psi)}.
  }
  \item{ACE.ci}{
    array; confidence interval determined by \code{\link{quantile}} if
    \code{ci.method} includes \dQuote{bootstrap}.  Otherwise calculated
    using analytic variance with large sample normal approximation.
    Array dimensions the same as \code{ACE} element.
  }
  \item{ACE.var}{
    array; estimated variance of ACE. Array dimensions the same as \code{ACE}
    element.
  }
  \item{ACE.p}{
    vector; estimated p-value of ACE.
  }
  \item{beta0}{
    vector; \eqn{\beta} values used for the first group.
  }
  \item{alphahat0}{
    vector; estimated \eqn{\alpha} values for the first group.
  }
  \item{Fas0}{
    function; estimator for the distribution function of \var{y0} in
    the first group in the always selected stratum.
  }
  \item{beta1}{
    vector; \eqn{\beta} values used for the second group.
  }
  \item{alphahat1}{
    vector; estimated \eqn{\alpha} values for the second group.
  }
  \item{Fas1}{
    function; estimator for the distribution function of \var{y1} in
    the second group in the always selected stratum.
  }
  \item{phi}{
    vector; \eqn{phi} values used.
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
  Jemiai Y (2005), \dQuote{Semiparametric Methods for Inferring Treatment
  Effects on Outcomes Defined Only if a Post-Randomization Event
  Occurs,} unpublished doctoral dissertation under the supervision of
  A. Rotnitzky, Harvard School of Public Health, Dept. of Biostatistics.

  Shepherd BE, Redman MW, Ankerst DP (2008), \dQuote{Does Finasteride affect the
  severity of prostate cancer? A causal sensitivity analysis,} Journal
  of the American Statistical Association 2008, 484, 1392-1404.
  
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
  \code{\link[sensitivityPStrat]{sensitivityGBH}},
  \code{\link[sensitivityPStrat]{sensitivitySGD}}
}
\examples{


data(vaccine.trial)
ansJR<-with(vaccine.trial,
          sensitivityJR(z=treatment,s=hiv.outcome,y=logVL,
                    beta0=c(-1,-.5,0,.5,1),
                    beta1=c(-1,-.5,0,.5,1),
                    phi=c(0.95,0.9), selection="infected",
                    groupings=c("placebo","vaccine"),
                    N.boot=100)
         )
ansJR

data(vaccine.trial)
ansJR<-with(vaccine.trial,
          sensitivityJR(z=treatment,s=hiv.outcome,y=logVL,
                    beta0=c(-1,-.5,0,.5,1),
                    beta1=c(-1,-.5,0,.5,1),
                    phi=c(0.95,0.9), selection="infected",
                    groupings=c("placebo","vaccine"),
                    custom.FUN=function(mu0, mu1, ...) mu1 - mu0,
                    upperTest=TRUE, lowerTest=TRUE, twoSidedTest=TRUE,
                    N.boot=100)
         )
ansJR

}

