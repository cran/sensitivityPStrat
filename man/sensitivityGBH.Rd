\name{sensitivityGBH}
\alias{sensitivityGBH}
\title{
  Principal stratification sensitivity analysis.
}
\description{
  Performs a sensitivity analysis using the method described in Gilbert,
  Bosch, and Hudgens (2003).
}
\usage{
sensitivityGBH(z, s, y, beta, selection, groupings,
               empty.principal.stratum, ci = 0.95,
               ci.method = c("analytic", "bootstrap"),
               ci.type = "twoSided", custom.FUN = NULL, na.rm = FALSE,
               N.boot = 100, interval = c(-100, 100),
               upperTest = FALSE, lowerTest = FALSE, twoSidedTest = TRUE,
               method = c("ACE", "T1", "T2"), isSlaveMode=FALSE)
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
    vector; outcome value. Can be \code{NA} for unselected records.
  }
  \item{beta}{
    vector; values of the \eqn{\beta} sensitivity parameter.  \code{Inf} and
    \code{-Inf} are acceptable.
  }
  \item{selection}{
    The value of \code{s} indicating selection.
  }
  \item{groupings}{
    vector of two elements \code{c(\var{g0},\var{g1})}; describes the
    possible group values. The first element \var{g0} being the value of
    \code{z} that delineates the first group, the last element \var{g1}
    being the value of \code{z} that delineates the second group.
  }
  \item{empty.principal.stratum}{
    vector of two elements \code{c(\var{s0},\var{s1})}; describes the \code{s}
    values that select the empty principal stratum.  If
    \code{empty.principal.stratum=c(\var{s0},\var{s1})}, then stratum defined by
    \eqn{S(\var{g0}) = \var{s0}} and \eqn{S(\var{g1}) = \var{s1}} is the empty
    stratum.  In this example \var{s0} and \var{s1} refer to the two
    possible values of \code{s}. (Note: method only works if \eqn{\var{s0}
    \ne \var{s1}}{\var{s0} != \var{s1}}).
  }
  \item{ci}{
    numeric vector; confidence interval level. Defaults to \code{0.95}
  }
  \item{ci.method}{
    character; method by which the confidence interval and variance are
    calculated.  Can be \dQuote{analytic} or
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
    \eqn{\code{p1} = P(S(\var{g1}) = selected)}. The custom function
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
    by \code{\link{optimize}} to estimate \eqn{\alpha}.
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
    character vector; type of test statistic calculated. Can be one or more
    of \dQuote{ACE}, \dQuote{T1}, or \dQuote{T2}. Defaults to
    \code{"ACE"}.  Methods \dQuote{T1} and \dQuote{T2} are not implemented
    if \code{ci.method} includes \dQuote{analytic}.
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
  empty).  ACE is identified by assuming a value of the sensitivity
  parameter beta, where \eqn{e^{\beta}}{exp(\beta)} has an odds ratio interpretation:

  If \eqn{\code{empty.principal.stratum}=c(S(\var{g0})=not\ selected, S(\var{g1})=selected)} then
  given selected if assigned \var{g0}, the odds of being selected if assigned
  \var{g1} multiplicatively increase \eqn{e^{\beta}}{exp(\beta)} for
  every 1-unit increase in \eqn{Y(\var{g0})}.
  
  If \eqn{\code{empty.principal.stratum}=c(S(\var{g0})=selected,
  S(\var{g1})=not\ selected)} then given selected if assigned \var{g1}, the odds of being selected if assigned \var{g0} multiplicatively increase
  \eqn{e^{\beta}}{exp(\beta)} for every 1-unit increase in \eqn{Y(\var{g1})}.

  Specifying \code{beta}=\code{-Inf} or \code{beta}=\code{Inf} calls
  \code{\link{sensitivityHHS}}.

  T1 and T2 are rank-based analogs of ACE. See <REF TBD>.
}
\value{
  an object of class \code{sensitivity2d}.

  \item{ACE}{
    vector;
    \eqn{ACE = E(Y(\var{g1}) - Y(\var{g0})|S(\var{g1}) = S(\var{g0}) = \code{selection})}.
    Vector of the estimated ACE values for specified \code{beta}
    values. Only exists if \code{method} includes \dQuote{ACE}.
  }
  \item{ACE.ci}{
    array; confidence interval of ACE determined by
    quantiles of bootstrap if \code{ci.method} includes
    \dQuote{bootstrap}.  Otherwise calculated using analytic variance
    with large sample normal approximation. Only exists if
    \code{method} includes \dQuote{ACE}.
  }
  \item{ACE.var}{
    vector; estimated variance of ACE.  Only exists if
    \code{method} includes \dQuote{ACE}.
  }
  \item{ACE.p}{
    vector; estimated p-value of ACE.  Only exists if
    \code{method} includes \dQuote{ACE}.
  }
  \item{T1}{
    vector; Vector of the estimated T1 test statistic for specified
    \code{beta} values. Only exists if \code{method} includes
    \dQuote{T1}.
  }
  \item{T1.p}{
    vector; estimated p-value of T1.  Only exists if
    \code{method} includes \dQuote{T1}.
  }
  \item{T2}{
    vector; Vector of the estimated T2 statistic for specified
    \code{beta} values. Only exists if \code{method} includes
    \dQuote{T2}.
  }
  \item{T2.p}{
    vector; estimated p-value of T2.  Only exists if
    \code{method} includes \dQuote{T2}.
  }
  \item{beta}{
    vector; user-specified \eqn{\beta} values
  }
  \item{alphahat}{
    vector; estimated values of \eqn{\alpha}
  }
  \item{Fas0}{
    function; estimator for the empirical distribution function values
    for \var{y0} in the first group in the always selected principal
    stratum. \eqn{Pr(Y(\var{g0}) \le
    \var{y0}|S(\var{g0})=S(\var{g1})=\code{selection};
    \beta)}{Pr(Y(\var{g0}) <= \var{y0}|S(\var{g1}) = \code{selection}; \beta)}
  }
  \item{Fas1}{
    function; estimator for the empirical distribution function values
    for \var{y1} in the second group in the always selected principal
    stratum.  \eqn{Pr(Y(\var{g1}) \le \var{y1}|S(\var{g0})=S(\var{g1})=\code{selection}; \beta)}{Pr(Y(\var{g1}) <= \var{y1}|S(\var{g0}) = S(\var{g1}) = \code{selection}; \beta)}
  }
}
\references{
  Gilbert PB, Bosch RJ, and Hudgens MG (2003), \dQuote{Sensitivity Analysis for
  the Assessment of Causal Vaccine Effects of Viral Load in HIV Vaccine
  Trials,} Biometrics 59, 531-541.
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
  \code{\link{sensitivityHHS}}, \code{\link{sensitivityJR}}, \code{\link{sensitivitySGL}}
}
\examples{
data(vaccine.trial)
ans<-with(vaccine.trial,
          sensitivityGBH(z=treatment,s=hiv.outcome,y=logVL,
                    beta=c(0,.25,.5,.75,1,1.25,1.5),
                    selection="infected",
                    groupings=c("placebo","vaccine"),
                    empty.principal.stratum=c("not infected","infected"),
                    N.boot=100)
         )
ans

ans<-with(vaccine.trial,
          sensitivityGBH(z=treatment,s=hiv.outcome,y=logVL,
                    beta=c(-Inf,-1,-0.75,-0.5,-0.25,0,.25,.5,.75,1,Inf),
                    selection="infected",
                    groupings=c("placebo","vaccine"),
                    empty.principal.stratum=c("not infected","infected"),
                    ci.method="bootstrap", ci=c(0.95, 0.9, 0.9),
                    ci.type=c('twoSided', 'upper', 'lower'),
                    custom.FUN=function(mu0, mu1, ...) mu1 - mu0,
                    N.boot=100, method=c("ACE", "T1", "T2"),
                    upperTest=TRUE, lowerTest=TRUE, twoSidedTest=TRUE)
         )
ans
}
\keyword{design}

