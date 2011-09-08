\name{calc.v}
\alias{calc.v}
\title{
  Calculates the v matrix used in the estimation of standard errors in sensitivitySGL.
}
\description{
  Calculates the v matrix used in the estimation of standard errors in sensitivitySGL.
}
\usage{
calc.v(event, time)
}
\arguments{
  \item{event}{
    logical vector indicating whether and event has happened.
  }
  \item{time}{
    vector; time until event or observation halted.
  }
}
\value{
  returns a matrix.
}
\references{
  Shepherd BE, Gilbert PB, Lumley T (2007), \dQuote{Sensitivity analyses
  comparing time-to-event outcomes existing only in a subset selected
  postrandomization,} Journal of the American Statistical Association
  102, 573-582.
}
\author{
  Bryan E. Shepherd \cr
  Department of Biostatistics\cr
  Vanderbilt University\cr

  
  Charles Dupont \cr
  Department of Biostatistics\cr
  Vanderbilt University\cr
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{design}
