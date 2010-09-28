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
  A treatment effect is a contrast between Y(0) and Y(1) where Y(0) is
  the outcome if not treated and Y(1) is the outcome if treated.  The
  average treatment effect (or average causal effect) is E(Y(1)-Y(0)).
  In some settings there may be interest in estimating the average
  treatment effect among those who would be selected under either
  treatment assignment (i.e., E(Y(1)-Y(0)|S(0)=S(1)=1), where S(0) is
  the indicator of selection if not treated and S(1) is the indicator of
  selection if treated (Robins 1986).  For example, one may want to
  assess the average treatment effect of a drug on quality of life among
  those who would have lived regardless of their treatment assignment.
  The subgroup defined by S(0)=S(1)=1 (e.g., those who would have lived
  regardless of treatment assignment) has been referred to as a
  principal stratum (Frangakis and Rubin, 2002).  Principal stratum
  membership is not known so to identify the average treatment effect
  (or related estimands) within a principal stratum we assume 1. SUTVA
  (Rubin 1978) (i.e., no interference -- that the potential outcomes for
  all subjects are independent of the treatment assignment of other
  subjects),2. ignorable treatment assignment (i.e., random assignment
  of treatment), 3. that one of the principal strata is empty, and
  4. that outcomes if assigned one treatment and selected are
  independent of selection if assigned the other treatment.  This
  package implements sensitivity analysis methods that relax these
  latter two assumptions.

  sensitivityHHS and sensitivityGBH implement the methods described by Hudgens, Hoering and Self (2003) and Gilbert, Bosch, and Hudgens (2003), respectively.  They estimate the average treatment effect in the always-selected principal stratum under assumptions 1-3, relaxing 4 using a worse-case scenario analysis (sensitivityHHS) or using a sensitivity parameter (sensitivityGBH).
  sensitivitySGL implements the methods described by Shepherd, Gilbert,
  and Lumley (2006).  It is similar to sensitivityHHS and sensitivityGBH
  except that it computes the difference between distribution functions
  in the always-selected principal stratum and allows the outcome to be
  right-censored.  sensitivityJR estimates the average treatment effect
  in the always-selected principal stratum relaxing assumptions 3 and 4
  as described by Jemiai and Rotnitzky (2005) and Shepherd, Redman, and
  Ankerst (2008).  sensitivitySGD incorporates the methods of Shepherd,
  Gilbert, and Dupont (in press), extending sensitivityJR to right-censored outcomes.
  
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
  Frangakis CE and Rubin DB (2002), "Principal stratification in causal inference," Biometrics 58, 21-29.

  Gilbert PB, Bosch RJ, and Hudgens MG (2003), "Sensitivity Analysis for
  the Assessment of Causal Vaccine Effects of Viral Load in HIV Vaccine
  Trials," Biometrics 59, 531-541.
  
  Hudgens MG, Hoering A, and Self SG (2003), "On the Analysis of Viral
  Load Endpoints in HIV Vaccien Trials," Statistics in Medicine 22,
  2281-2298.

  Jemiai Y (2005), "Semiparametric Methods for Inferring Treatment
  Effects on Outcomes Defined Only if a Post-Randomization Event
  Occurs," unpublished doctoral dissertation under the supervision of
  A. Rotnitzky, Harvard School of Public Health, Dept. of Biostatistics.
  
  Robins JM (1986), "A new approach to causal inference in mortality studies with sustained exposure periods - Application to control of the healthy worker survivor effect," Mathemeatical Modeling 7, 1393-1512.

  Rubin DB (1978), "Bayesian inference for causal effects: the role of randomization," The Annals of Statistics 6, 34-58.

  Shepherd BE, Gilbert PB, Lumley T (2007), "Sensitivity analyses
  comparing time-to-event outcomes existing only in a subset selected
  postrandomization," Journal of the American Statistical Association
  102, 573-582.
  
  Shepherd BE, Gilbert PB, and Dupont CT, "Sensitivity analyses
  comparing time-to-event outcomes only existing in a subset selected
  postrandomization and relaxing monotonicity," Biometrics (in press).
}
\keyword{ package }
\seealso{
  \code{\link[survival]{Surv}}
}
