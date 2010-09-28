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
              na.rm = FALSE, N.boot = 100, oneSidedTest = FALSE,
              twoSidedTest = TRUE, verbose=getOption("verbose"))
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
    outcome vector.  Can be \code{NA} for
    unselected records.
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
  \item{selection}{
    The value of \code{s} indicating selection.
  }
  \item{groupings}{
    Vector of two elements \code{c(g0,g1)}, the first element \code{g0} being the
    value of \code{z} the delineates the first group, the last element
    \code{g1} being the value of \code{z} which delineates the second group.
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
  \item{verbose}{
    logical;  prints dots when bootstrapping to show that something is happening.
  }
}
\details{
  
  Performs a sensitivity analysis estimating the average causal effect
  among those who would have been selected regardless of treatment
  assignment (ACE) without assuming monotonicity (i.e., that one of the
  principal strata is empty).  The method assumes no interference (i.e., potential outcomes of all
  subjects are unaffected by treatment assignment of other subjects) and
  ignorable (i.e., random) treatment
  assignment.  ACE is identified by assuming values for the sensitivity parameters beta0,
    beta1, and one of the parameters phi, psi, or Pi.  The sensitivity
    parameters beta0 and beta1 have a log-odds ratio interpretation (see
    help for sensitivityGBH).  
  
  Only one of the parameters phi, psi, or Pi should be specified as all
  depend on each other.  psi is unrestrained taking any value on the real
  line.  The other parameters, psi and Pi have constraints and there
  will be estimation problems if these parameters are set at values
  outside the of their range of acceptable values based on the observed
  data.  See Shepherd, Gilbert, Dupont (in press) for more details.

}
\value{
  object of class \code{sensitivity3d}

  \item{ACE}{
    array; estimated values of ACE for all combinations of \code{beta0},
    \code{beta1}, and \code{phi/Pi/psi}.  Array dimensions are \code{length(beta0)},
    \code{length(beta1)}, \code{length(psi)}.
  }
  \item{beta0}{
    beta0 values used.
  }
  \item{beta1}{
    beta1 values used.
  }
  \item{phi}{
    phi values used.
  }
  \item{ACE.ci}{
    array; confidence interval determined by \code{quantile} if using
    \code{ci.method} \dQuote{bootstrap}.  Otherwise calculated using
    analytic variance with large sample normal approximation.  Array
    dimensions the same as \code{ACE} element.
  }
  \item{ACE.var}{
    array; estimated variance of ACE.   Array dimensions the same as \code{ACE}
    element.
  }
}
\references{
  Jemiai Y (2005), "Semiparametric Methods for Inferring Treatment
  Effects on Outcomes Defined Only if a Post-Randomization Event
  Occurs," unpublished doctoral dissertation under the supervision of
  A. Rotnitzky, Harvard School of Public Health, Dept. of Biostatistics.

  Shepherd BE, Redman MW, Ankerst DP (2008), "Does Finasteride affect the
  severity of prostate cancer? A causal sensitivity analysis," Journal
  of the American Statistical Association 2008, 484, 1392-1404.
  
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
  \code{\link[sensitivityPStrat]{sensitivityGBH}},
  \code{\link[sensitivityPStrat]{sensitivitySGD}}
}
\examples{


data(vaccine.trial)
ans<-with(vaccine.trial,
          sensitivityJR(z=treatment,s=hiv.outcome,y=logVL,
                    beta0=c(0,.25,.5,.75,1), beta1=c(0,.25,.5,.75,1),
                    phi=c(0.95,0.90,0.80), selection="infected",
                    groupings=c("placebo","vaccine"),
                    N.boot=1000)
         )
ans

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

