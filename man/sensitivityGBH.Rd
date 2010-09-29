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
               ci.method = c("analytic", "bootstrap"), na.rm = FALSE,
               N.boot = 100, oneSidedTest = FALSE, twoSidedTest = TRUE)
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
  \item{beta}{
    vector; values of the sensitivity parameter.  \code{Inf} and
    \code{-Inf} are acceptable.
  }
  \item{selection}{
    The value of \code{s} indicating selection.
  }
  \item{groupings}{
    Vector of two elements \code{c(g0,g1)}, the first element \code{g0} being the
    value of \code{z} the delineates the first group, the last element
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
    \dQuote{bootstrap}. Defaults to \code{c("analytic","bootstrap")}
  }
  \item{na.rm}{
    logical; indicates whether records that are invalid due to \code{NA}
    values should be
    removed from the data set.
  }
  \item{N.boot}{
    integer; number of bootstrap repetitions that will be run when
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
}
\details{
  Performs a sensitivity analysis estimating the average causal effect
  among those who would have been selected regardless of treatment
  assignment (ACE).  The method assumes no interference (i.e., potential outcomes of all
  subjects are unaffected by treatment assignment of other subjects),
  ignorable (i.e., random) treatment
  assignment, and monotonicity (i.e., one of the principal strata is
  empty).  ACE is identified by assuming a value of the sensitivity
  parameter beta, where exp(beta) has an odds ratio interpretation:

  If empty.principal.stratum=c(S(g0)==not selected, S(g1)==selected) then
  given selected if assigned g0, the odds of being selected if assigned
  g1 multiplicatively increase exp(beta) for every 1-unit increase in
  Y(g0).

  If empty.principal.stratum=c(S(g0)==selected, S(g1)==not
  selected) then given selected if assigned g1, the odds of being
  selected if assigned g0 multiplicatively increase exp(beta) for every
  1-unit increase in Y(g1).

  Specifying beta=\code{-Inf} or
  beta=\code{Inf} calls \code{\link{sensitivityHHS}}.
  }
\value{
  an object of class \code{sensitivity2d}.

  \item{ACE}{
    ACE=E(Y(g1)-Y(g0)|S(g1)==S(g0)==selection).  Vector of the estimated ACE values for specified \code{beta} values.
  }
  \item{ACE.ci}{
    vector; confidence interval of ACE determined by
    quantiles of bootstrap if \code{ci.method} includes
    \dQuote{bootstrap}.  Otherwise calculated using analytic variance
    with large sample normal approximation.
  }
  \item{ACE.var}{
    vector; estimated variance of ACE.
  }
  \item{beta}{
    vector of user-specified beta values
  }
  \item{alphahat}{
    vector of estimated values of alpha
  }
  \item{y0}{
    vector of unique \code{y} values in the first group.
  }
  \item{Fas0}{
    matrix of estimated empirical distribution function values for
    \code{y0} in the first group in the always selected principal
    stratum. Pr(Y(g0) <= y0|S(g0)=S(g1)=selection; beta)
  }
  \item{y1}{
    vector of unique \code{y} values in the second group.
  }
  \item{Fas1}{
    matrix of estimated empirical distribution function values for
    \code{y1} in the second group in the always selected principal
    stratum.  Pr(Y(g1) <= y1|S(g0)=S(g1)=selection; beta)
  }
}
\references{
  Gilbert PB, Bosch RJ, and Hudgens MG (2003), "Sensitivity Analysis for
  the Assessment of Causal Vaccine Effects of Viral Load in HIV Vaccine
  Trials," Biometrics 59, 531-541.
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
                    N.boot=1000)
         )
ans
\dontshow{
stopifnot(is.list(ans))
stopifnot(inherits(ans,"sensitivity2d"))
}

ans<-with(vaccine.trial,
          sensitivityGBH(z=treatment,s=hiv.outcome,y=logVL,
                    beta=c(-Inf, 0,.25,.5,.75,1,Inf),                      
                    selection="infected",
                    groupings=c("placebo","vaccine"),
                    empty.principal.stratum=c("not infected","infected"),
                    ci.method="bootstrap",
                    N.boot=1000)
         )
ans
\dontshow{
stopifnot(is.list(ans))
stopifnot(inherits(ans,"sensitivity2d"))
}
ans<-with(vaccine.trial,
          sensitivityGBH(z=treatment,s=hiv.outcome,y=logVL,
                    beta=-Inf,                      
                    selection="infected",
                    groupings=c("placebo","vaccine"),
                    empty.principal.stratum=c("not infected","infected"),
                    ci.method="bootstrap",
                    N.boot=1000)
         )
ans
\dontshow{
stopifnot(is.list(ans))
stopifnot(inherits(ans,"sensitivity2d"))
}
}
\keyword{design}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

