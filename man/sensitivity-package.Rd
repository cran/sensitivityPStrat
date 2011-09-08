\name{sensitivityPStrat-package}
\alias{sensitivityPStrat-package}
\alias{sensitivityPStrat}
\docType{package}
\title{
  Principal Stratification Sensitivity Analysis Functions
}
\description{
  This package provides functions to perform sensitivity analyses of treatment effects within principal strata.
}
\details{
  A treatment effect is a contrast between \eqn{Y(0)} and \eqn{Y(1)} where \eqn{Y(0)} is
  the outcome if not treated and \eqn{Y(1)} is the outcome if treated.  The
  average treatment effect (or average causal effect) is \eqn{E(Y(1)-Y(0))}.
  In some settings there may be interest in estimating the average
  treatment effect among those who would be selected under either
  treatment assignment (i.e., \eqn{E(Y(1)-Y(0)|S(0)=S(1)=1)}, where \eqn{S(0)} is
  the indicator of selection if not treated and \eqn{S(1)} is the indicator of
  selection if treated (Robins 1986).  For example, one may want to
  assess the average treatment effect of a drug on quality of life among
  those who would have lived regardless of their treatment assignment.
  The subgroup defined by \eqn{S(0)=S(1)=1} (e.g., those who would have lived
  regardless of treatment assignment) has been referred to as a
  principal stratum (Frangakis and Rubin, 2002).  Principal stratum
  membership is not known so to identify the average treatment effect
  (or related estimands) within a principal stratum we assume 1. SUTVA
  (Rubin 1978) (i.e., no interference -- that the potential outcomes for
  all subjects are independent of the treatment assignment of other
  subjects),2. ignorable treatment assignment (i.e., random assignment
  of treatment), 3. that one of the principal strata is empty, and
  4. that a selected subject's outcome if assigned one treatment is
  independent of selection if assigned the other treatment.  This
  package implements sensitivity analysis methods that relax these
  latter two assumptions.

  \code{\link{sensitivityHHS}} and \code{\link{sensitivityGBH}}
  implement the methods described by Hudgens, Hoering and Self (2003)
  and Gilbert, Bosch, and Hudgens (2003), respectively.  They estimate
  the average treatment effect in the always-selected principal stratum
  under assumptions 1-3, relaxing 4 using a worse-case scenario analysis
  (\code{sensitivityHHS}) or using a sensitivity parameter
  (\code{sensitivityGBH}).  These functions also have options to do
  rank-based analyses and to compute other measures of treatment
  efficacy with continuous or binary outcomes (Hudgens and Halloran, 2006).
  \code{\link{sensitivitySGL}} implements the methods described by Shepherd, Gilbert,
  and Lumley (2006).  It is similar to \code{sensitivityHHS} and \code{sensitivityGBH}
  except that it computes the difference between distribution functions
  in the always-selected principal stratum and allows the outcome to be
  right-censored.  \code{\link{sensitivityJR}} estimates the average treatment effect
  in the always-selected principal stratum relaxing assumptions 3 and 4
  as described by Jemiai and Rotnitzky (2005) and Shepherd, Redman, and
  Ankerst (2008).  \code{\link{sensitivitySGD}} incorporates the methods of Shepherd,
  Gilbert, and Dupont (in press), extending \code{sensitivityJR} to right-censored outcomes.
  
}
\author{
  Bryan E. Shepherd \cr
  Department of Biostatistics\cr
  Vanderbilt University\cr

  
  Charles Dupont \cr
  Department of Biostatistics\cr
  Vanderbilt University\cr
  \email{charles.dupont@vanderbilt.edu}

  Maintainer: Charles Dupont
}
\references{
  Frangakis CE and Rubin DB (2002), \dQuote{Principal stratification in causal inference,} Biometrics 58, 21-29.

  Gilbert PB, Bosch RJ, and Hudgens MG (2003), \dQuote{Sensitivity Analysis for
  the Assessment of Causal Vaccine Effects of Viral Load in HIV Vaccine
  Trials,} Biometrics 59, 531-541.

  Hudgens MG, Halloran ME, \dQuote{Causal vaccine effects on binary post
  infection outcomes,} Journal of the American Statisitcal Association
  101, 51-64. 
  
  Hudgens MG, Hoering A, and Self SG (2003), \dQuote{On the Analysis of Viral
  Load Endpoints in HIV Vaccine Trials,} Statistics in Medicine 22,
  2281-2298.

  Jemiai Y (2005), \dQuote{Semiparametric Methods for Inferring Treatment
  Effects on Outcomes Defined Only if a Post-Randomization Event
  Occurs,} unpublished doctoral dissertation under the supervision of
  A. Rotnitzky, Harvard School of Public Health, Dept. of Biostatistics.
  
  Robins JM (1986), \dQuote{A new approach to causal inference in mortality studies with sustained exposure periods - Application to control of the healthy worker survivor effect,} Mathematical Modeling 7, 1393-1512.

  Rubin DB (1978), \dQuote{Bayesian inference for causal effects: the role of randomization,} The Annals of Statistics 6, 34-58.

  Shepherd BE, Gilbert PB, Lumley T (2007), \dQuote{Sensitivity analyses
  comparing time-to-event outcomes existing only in a subset selected
  postrandomization,} Journal of the American Statistical Association
  102, 573-582.
  
  Shepherd BE, Gilbert PB, and Dupont CT, \dQuote{Sensitivity analyses
  comparing time-to-event outcomes only existing in a subset selected
  postrandomization and relaxing monotonicity,} Biometrics (in press).
}
\keyword{ package }
\seealso{
  \code{\link[survival]{Surv}}
}
