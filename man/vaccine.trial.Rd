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
    \item{\code{treatment}}{a factor with levels \code{placebo} \code{vaccine}}
    \item{\code{hiv.outcome}}{a factor with levels \code{infected} \code{not infected}}
    \item{\code{logVL}}{a numeric vector}
    \item{\code{ARTinitiation}}{a factor with levels \code{initiated ART} \code{no ART}}
    \item{\code{followup.yearsART}}{a numeric vector}
  }
}
\examples{
\dontrun{
# Put code which generates the simluated dataset here
}
}
\keyword{datasets}
