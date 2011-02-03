\name{sensitivityHHS}
\alias{sensitivityHHS}
\title{
  principal stratifictation sensitivity analysis using the HHS method.
}
\description{
  Performs a principal stratifictation sensitivity analysis using the
  method described in Hudgens, Hoering, and Self (2003).
}
\usage{
sensitivityHHS(z, s, y, bound = c("upper", "lower"), selection,
               groupings, empty.principal.stratum, ci = 0.95,
               ci.method = c("bootstrap", "analytic"), na.rm = FALSE,
               N.boot = 100, oneSidedTest = FALSE, twoSidedTest = TRUE,
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
    outcome vector.  Can be \code{NA} for
    unselected records.
  }
  \item{bound}{
    vector selecting which bound should be calculated, upper and/or lower.  Partial string
    matching is performed.
  }
  \item{selection}{
    The value of \code{s} indicating selection.
  }
  \item{groupings}{
    Vector of two elements \code{c(g0,g1)}, first element \code{g0} being the
    value of \code{z} which delineates the first group, the last element
    \code{g1} being the value of \code{z} which delineates the second group.
  }
  \item{empty.principal.stratum}{
    vector of two elements \code{c(s0,s1)}; describes the \code{s}
    values that select the empty principal stratum.  If
    \code{empty.principal.stratum=c(s0,s1)}, then stratum defined by
    S(g0)==s0 and S(g1)==s1 is the empty stratum.  In this example s0
    and s1 refer to the two possible values of s. (Note: method only
    works if s0 != s1).
  }
  \item{ci}{
    numeric vector;  confidence interval level, defaults to 0.95
  }
  \item{ci.method}{
    character;  method by which the confidence interval and
    variance are calculated.  Can be \dQuote{analytic} or
    \dQuote{bootstrap}. Defaults to \code{c("analytic","bootstrap")}.
    Currently only works for \dQuote{bootstrap}.
  }
  \item{na.rm}{
    logical; indicates whether records that are invalid due to \code{NA}
    values should be
    removed from the data set.
  }
  \item{N.boot}{
    integer. Number of bootstrap repetitions that will be run when
    \code{ci.method} includes \dQuote{bootstrap}.
  }
  \item{oneSidedTest}{
    logical.  Return a one sided confidence interval for ACE. Defaults
    to \code{FALSE}
  }
  \item{twoSidedTest}{
    logical.  Return a two sided confidence interval for ACE. Defaults
    to \code{TRUE}
  }
  \item{isSlaveMode}{
    logical.  Internal Use only. Used in recursion.
  }
}
\details{
  Performs a sensitivity analysis estimating the average causal effect
  among those who would have been selected regardless of treatment
  assignment (ACE).  The method assumes no interference (i.e., potential outcomes of all
  subjects are unaffected by treatment assignment of other subjects),
  ignorable (i.e., random) treatment
  assignment, and monotonicity (i.e., one of the principal strata is
  empty).  ACE is still not identified after making these assumptions,
  so this method computes the lower and upper bounds of the estimated
  ACE.  These bounds correspond to the values one would get if using
  \code{\link{sensitivityGBH}} and specifying the sensitivity parameter beta as
  \code{-Inf} or \code{Inf}.
  }
\value{
  an object of class \code{sensitivity2d}.
  
  \item{ACE}{
    ACE=E(Y(g1)-Y(g0)|S(g1)==S(g0)==selection).  Vector of the estimated
    ACE values at the specified bounds.
  }
  \item{ACE.ci}{
    vector; confidence interval of ACE determined by
    quantiles of bootstrap if \code{ci.method} includes
    \dQuote{bootstrap}.  Otherwise calculated using analytic variance
    with large sample normal approximation (NOT YET WORKING).
  }
  \item{ACE.var}{
    vector; estimated variance of ACE.
  }
  \item{y0}{
    vector of unique \code{y} values in the first group.
  }
  \item{Fas0}{
    matrix of estimated empirical distribution function values for
    \code{y0} in the first group in the always selected principal
    stratum at the bounds. Pr(Y(g0) <= y0|S(g0)=S(g1)=selection)
  }
  \item{y1}{
    vector of unique \code{y} values in the second group.
  }
  \item{Fas1}{
    matrix of estimated empirical distribution function values for
    \code{y1} in the second group in the always selected principal
    stratum at the bounds.  Pr(Y(g1) <= y1|S(g0)=S(g1)=selection)
  }
}
\references{
  Hudgens MG, Hoering A, and Self SG (2003), "On the Analysis of Viral
  Load Endpoints in HIV Vaccien Trials," Statistics in Medicine 22,
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
}

