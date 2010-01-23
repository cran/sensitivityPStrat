funArray <- function(data=NA, dim=length(data), dimnames=NULL) {
  fcall <- match.call()

  fcall[[1]] <- as.name('array')
  
  data <- eval.parent(fcall)
  
  class(data) <- 'funArray'

  data
}

funVector <- function(length = 0) {
  x <- vector('list', length)
  class(x) <- 'funVector'

  x
}

'[.funVector' <- function(x, ..., drop=TRUE) {
  cl <- oldClass(x)
  class(x) <- NULL
  
  x <- NextMethod(.Generic)

  if(length(x) == 1)
    return(x[[1]])
 
  class(x) <- cl
  x
}

'[.funArray' <- function(x, ..., drop=TRUE) {
  cl <- oldClass(x)
  class(x) <- NULL
  
  x <- NextMethod(.Generic)

  if(is.array(x))
    class(x) <- cl
  else if(length(x) == 1)
    return(x[[1]])
  else
    class(x) <- ifelse(cl == 'funArray', 'funVector', cl)
  
  x
}

'[<-.funVector' <- function(x, ..., value) {
  if (!as.logical(length(value)))
    return(x)

  cl <- oldClass(x)
  if(is.function(value))
    value <- list(value)
  
  class(x) <- class(value) <- NULL

  x <- NextMethod(.Generic)

  class(x) <- cl
  x
}

'[<-.funArray' <- function(x, ..., value) {
  if (!as.logical(length(value)))
    return(x)
  
  cl <- oldClass(x)
  if(is.function(value))
    value <- list(value)
  
  class(x) <- class(value) <- NULL

  x <- NextMethod(.Generic)

  class(x) <- cl
  x
}

