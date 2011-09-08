\name{vaccine.trial}
\alias{vaccine.trial}
\docType{data}
\title{
  Simulated Vaccine Trial Data
}
\description{
  Simulated vaccine trial data for use in demonstrating the use of the
  sensitivity functions implemented in this package.
}
\usage{data(vaccine.trial)}
\format{
  A data frame with 2000 observations on the following 5 variables.
  \describe{
    \item{\code{treatment}}{a factor with levels \dQuote{placebo} ,\dQuote{vaccine}}
    \item{\code{hiv.outcome}}{a factor with levels \dQuote{infected}, \dQuote{not infected}}
    \item{\code{logVL}}{a numeric vector}
    \item{\code{ARTinitiation}}{a factor with levels \dQuote{initiated ART}, \dQuote{no ART}}
    \item{\code{followup.yearsART}}{a numeric vector}
  }
}
\keyword{datasets}
