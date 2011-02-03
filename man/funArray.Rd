\name{funArray}
\alias{funArray}
\title{
  Create an array of functions
}
\description{
  Creates a array of functions.
}
\usage{
funArray(data = NA, dim = length(data), dimnames = NULL)
}
\arguments{
  \item{data}{
    funVector; vector to turn into an array.
  }
  \item{dim}{
    integer vector; vector of maximal indices for each dimension.
  }
  \item{dimnames}{
    list of character vectors; either \code{NULL} or the names for each
    dimension. One component for each dimension, either \code{NULL} or
    the names for each dimension index.
  }
}
\author{
  Charles Dupont \cr
  Department of Biostatistics\cr
  Vanderbilt University\cr  
}
\seealso{
  \code{\link{funVector}}, \code{\link{vector}}
}
\keyword{ classes }
