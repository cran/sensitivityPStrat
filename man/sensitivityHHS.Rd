\name{sensitivityHHS}
\alias{sensitivityHHS}
\title{
  principal stratification sensitivity analysis using the HHS method.
}
\description{
  Performs a principal stratification sensitivity analysis using the
  method described in Hudgens, Hoering, and Self (2003).
}
\usage{
sensitivityHHS(z, s, y, bound = c("upper", "lower"), selection,
               groupings, empty.principal.stratum, ci = 0.95,
               ci.method = c("bootstrap", "analytic"),
               ci.type = "twoSided", custom.FUN = NULL, na.rm = FALSE,
               N.boot = 100, upperTest = FALSE, lowerTest = FALSE,
               twoSidedTest = TRUE, method = c("ACE", "T1", "T2"),
               isSlaveMode=FALSE)
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
  \item{y}{
    vector; outcome values. Can be \code{NA} for unselected records.
  }
  \item{bound}{
    vector; which bound should be calculated, \dQuote{upper} and/or
    \dQuote{lower}.  Partial string matching is performed.
  }
  \item{selection}{
    The value of \code{s} indicating selection.
  }
  \item{groupings}{
    vector of two elements \code{c(\var{g0},\var{g1})}; describes to
    possible group values. The first element \var{g0} being the value of
    \code{z} which delineates the first group, the last element \var{g1}
    being the value of \code{z} which delineates the second group.
  }
  \item{empty.principal.stratum}{
    vector of two elements \code{c(\var{s0},\var{s1})}; describes the
    \code{s} values that select the empty principal stratum.  If
    \code{empty.principal.stratum=c(\var{s0},\var{s1})}, then stratum
    defined by \eqn{S(\var{g0}) = \var{s0}} and
    \eqn{S(\var{g1}) = \var{s1}} is the empty stratum.  In this example
    \var{s0} and \var{s1} refer to the two possible values of
    \code{s}. (Note: method only works if
    \eqn{\var{s0} \ne \var{s1}}{\var{s0} != \var{s1}}).
  }
  \item{ci}{
    numeric vector;  confidence interval level, defaults to \code{0.95}.
  }
  \item{ci.method}{
    character;  method by which the confidence interval and
    variance are calculated.  Can be \dQuote{analytic} or
    \dQuote{bootstrap}. Defaults to \code{c("analytic","bootstrap")}.
    Currently only works for \dQuote{bootstrap}.
  }
  \item{ci.type}{
    character vector; type of confidence interval that the corisponding
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
    integer. Number of bootstrap repetitions that will be run when
    \code{ci.method} includes \dQuote{bootstrap}.
  }
  \item{lowerTest}{
    logical.  Return the lower one sided p-value for returned tests. Defaults
    to \code{FALSE}
  }
  \item{upperTest}{
    logical.  Return the upper one sided p-value for returned tests. Defaults
    to \code{FALSE}
  }
  \item{twoSidedTest}{
    logical.  Return a two sided p-value for returned tests. Defaults
    to \code{TRUE}
  }
  \item{method}{
    character vector; type of test statistic calculated. Can be one or
    more of \dQuote{ACE}, \dQuote{T1}, or \dQuote{T2}. Defaults to
    \code{"ACE"}.
  }
  \item{isSlaveMode}{
    logical; Internal Use only. Used in recursion.
  }
}
\details{
  Performs a sensitivity analysis estimating the average causal effect
  among those who would have been selected regardless of treatment
  assignment (ACE).  The method assumes no interference (i.e., potential outcomes of all subjects are unaffected by treatment assignment of other subjects),
  ignorable (i.e., random) treatment assignment, and monotonicity (i.e.,
  one of the principal strata is empty).  ACE is still not identified
  after making these assumptions, so this method computes the lower and
  upper bounds of the estimated ACE.  These bounds correspond to the
  values one would get if using \code{\link{sensitivityGBH}} and
  specifying the sensitivity parameter beta as \code{-Inf} or
  \code{Inf}.
}
\value{
  an object of class \code{sensitivity2d}.
  
  \item{ACE}{
    \eqn{ACE=E(Y(\var{g1})-Y(\var{g0})|S(\var{g1})=S(\var{g0})=\code{selection})}.
    Vector of the estimated ACE values at the specified bounds. Only
    exists if \code{method} includes \dQuote{ACE}.
  }
  \item{ACE.ci}{
    vector; confidence interval of ACE determined by
    quantiles of bootstrap if \code{ci.method} includes
    \dQuote{bootstrap}.  Otherwise calculated using analytic variance
    with large sample normal approximation (NOT YET WORKING). Only
    exists if \code{method} includes \dQuote{ACE}. 
  }
  \item{ACE.var}{
    vector; estimated variance of ACE. Only exists if \code{method}
    includes \dQuote{ACE}.
  }
  \item{ACE.p}{
    vector; estimated p-value of ACE.  Only exists if
    \code{method} includes \dQuote{ACE}.
  }
  \item{Fas0}{
    function; estimator for the empirical distribution function values
    for \var{y0} in the first group in the always selected principal
    stratum at the bounds.
    \eqn{Pr(Y(\var{g0}) \le \var{y0}|S(\var{g0})=S(\var{g1})=\code{selection})}{%
      Pr(Y(\var{g0}) <= \var{y0}|S(\var{g0}) = S(\var{g1}) = \code{selection})}
  }
  \item{Fas1}{
    function; estimator for the empirical distribution function values for
    \var{y1} in the second group in the always selected principal
    stratum at the bounds.
    \eqn{Pr(Y(\var{g1}) \le \var{y1}|S(\var{g0})=S(\var{g1})=\code{selection})}{%
      Pr(Y(\var{g1}) <= \var{y1}|S(\var{g0}) = S(\var{g1}) = \code{selection})}
  }
}
\references{
  Hudgens MG, Hoering A, and Self SG (2003), \dQuote{On the Analysis of Viral
  Load Endpoints in HIV Vaccine Trials,} Statistics in Medicine 22,
  2281-2298.
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
 \code{\link{sensitivityGBH}}, \code{\link{sensitivityJR}}, \code{\link{sensitivitySGL}}
}
\examples{
data(vaccine.trial)
est.bounds<-with(vaccine.trial,
                 sensitivityHHS(z=treatment, s=hiv.outcome, y=logVL,
                     selection="infected", groupings=c("placebo","vaccine"),
                     empty.principal.stratum=c("not infected","infected"),
                     N.boot=100)
                )
est.bounds

est.bounds<-with(vaccine.trial,
                 sensitivityHHS(z=treatment, s=hiv.outcome, y=logVL,
                     selection="infected", groupings=c("placebo","vaccine"),
                     empty.principal.stratum=c("not infected","infected"),
                     method=c("ACE", "T1", "T2"), N.boot=100,
                     custom.FUN=function(mu0, mu1, ...) mu1 - mu0,
                     upperTest=TRUE, lowerTest=TRUE, twoSidedTest=TRUE)
                )
est.bounds

}

