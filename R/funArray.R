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
  y <- NextMethod('[')

  cl <- oldClass(x)
  
  if(is.array(y))
    class(y) <- cl
  else if(length(y) == 1)
    return(y[[1]])
  else
    class(y) <- ifelse(cl == 'funArray', 'funVector', cl)
  
  y
}

## 'str.funVector' <- function(x, ...) {
##   P0 <- function(...) paste(..., sep = "")                                                                            
##   mod <- "function"

##   if (is.array(object)) {
##     rnk <- length(di. <- dim(object))
##     di <- P0(ifelse(di. > 1, "1:", ""), di., ifelse(di. >
##                                                     0, "", " "))
##     pDi <- function(...) paste(c("[", ..., "]"),
##                                collapse = "")
##     le.str <- (if (rnk == 1)
##                pDi(di[1L], "(1d)")
##     else pDi(P0(di[-rnk], ", "), di[rnk]))
##     std.attr <- "dim"
##   }
##   else if (!is.null(names(object))) {
##     mod <- paste("Named", mod)
##     std.attr <- std.attr[std.attr != "names"]
##   }
##   if (has.class && length(cl) == 1) {
##     if (cl != mod && substr(cl, 1, nchar(mod)) !=
##         mod)
##       mod <- P0("'", cl, "' ", mod)
##     std.attr <- c(std.attr, "class")
##   }
##   str1 <- if (le == 1 && !is.array(object))
##     paste(NULL, mod)
##   else P0(" ", mod, if (le > 0)
##           " ", le.str)


## }

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

