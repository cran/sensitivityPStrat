\name{sensitivityHHS}
\alias{sensitivityHHS}
\title{
  principle stratifictation sensitivity analysis using the HHS method.
}
\description{
  Preforms a principle stratifictation sensitivity analysis using the
  method described in HHS.
}
\usage{
sensitivityHHS(z, s, y, bound = c("upper", "lower"), selection,
               groupings, empty.principle.stratum, ci = 0.95,
               ci.method = c("bootstrap", "analytic"), na.rm = FALSE,
               N.boot = 100, oneSidedTest = FALSE, twoSidedTest = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
    vector containing the grouping values of the the records.
  }
  \item{s}{
    vector describing whether a record is selected.
  }
  \item{y}{
    vector describing the length after the event occured until observed
    effect.  Can be \code{NA} for unselected records.
  }
  \item{bound}{
    vector selecting which bound should be calculated.  Partial string
    matching is performed.
  }
  \item{selection}{
    The value of \code{s} in where a record is selected.
  }
  \item{groupings}{
    Vector of two elements \code{c(g0,g1)}, first element \code{g0} being the
    value of \code{z} which delineates the first group, the last element
    \code{g1} being the value of \code{z} which delineates the second group.
  }
  \item{empty.principle.stratum}{
    vector of two elements \code{c(s0,s1)}. Describing the \code{s}
    values which select the empty principle stratum.  If
    \code{empty.principle.stratum=c(s0,s1)}, then stratum defined by
    S(g0)==s0 and S(g1)==s1 is the empty stratum.  In this example s0
    and s1 refer to the two possible values of s. (Note: method only
    works if s0 != s1).
  }
  \item{ci}{
    numeric.  confindence interval level. defaults to 0.95
  }
  \item{ci.method}{
    character.  Method by which the convidence interval and
    variance are calculated.  Can be \dQuote{analytic} and
    \dQuote{bootstrap}. defaults to \code{c("analytic","bootstrap")}
  }
  \item{na.rm}{
    logical. Should records which are invalid due to \code{NA} values be
    removed from the data set.
  }
  \item{N.boot}{
    integer. Number of bootstrap repitions that will be run when
    \code{ci.method} includes \dQuote{bootstrap}.
  }
  \item{oneSidedTest}{
    logical.  Return a one sided confidence interval for ACE. defaults
    to \code{FALSE}
  }
  \item{twoSidedTest}{
    logical.  Return a two sided confidence interval for ACE. defaults
    to \code{TRUE}
  }
}
\value{
  an object of class \code{sensitivity2d}.
  
  \item{ACE}{
    vector of the calculated ACE values for the data set and
    bound selection
  }
  \item{ACE.ci}{
    vector; confidance interval of ACE determined by
    \code{\link{quantile}} if \code{ci.method} includes
    \dQuote{bootstrap}.  Other wise calculated using analytic varance
    and normal approximation.
  }
  \item{ACE.var}{
    vector; variance of ACE.
  }
  \item{beta}{
    beta values selected to ACE values.
  }
  \item{y0}{
    vector of unique \code{y} values in the first group.
  }
  \item{Fas0}{
    matrix of estimated emprical distribution function values for
    \code{y} in the first group in the always selected principle
    stratum. \code{S(g0) == selection}.
  }
  \item{y1}{
    vector of unique \code{y} values in the second group.
  }
  \item{Fas1}{
    matrix of estimated emprical distribution function values for
    \code{y} in the second group in the always selected principle
    stratum. \code{S(g1) == selection}.
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
  \code{\link[stats]{quantile}}, \code{\link[stats]{var}},
  \code{\link{qnorm}}, \code{\link{sensitivityGBH}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
\keyword{design}
