.CheckEmptyPrincipalStratum <- function(empty.principal.stratum, ...) {
  ## empty.principal.stratum must be
  ## 1. not missing
  ## 2. not NULL
  ## 3. a 2 element vector
  ## 4. elements must be diffrent values
  ## 5. cannot be NA

  if(missing(empty.principal.stratum))
    return("'empty.principal.stratum' argument must be supplied")

  ErrMsg <- NULL
  if(is.null(empty.principal.stratum) ||
     length(empty.principal.stratum) != 2L)
    ErrMsg <- c(ErrMsg,
                "'empty.principal.stratum' argument must be a two element vector")
  else if(empty.principal.stratum[1] == empty.principal.stratum[2])
    ErrMsg <- c(ErrMsg,
                "'empty.principal.stratum' argument cannot have the same values for both elements")
  

  if(length(empty.principal.stratum) &&
     any(is.na(empty.principal.stratum)))
    ErrMsg <- c(ErrMsg,
                "'empty.principal.stratum' may not contain a NA")

  return(ErrMsg)
}


.CheckSelection <- function(selection, s, empty.principal.stratum, ...) {
  ## selection must be
  ## 1. not missing
  ## 2. a single element vector
  ## 3. not NA
  ## 4. contained in empty.principal.stratum if it exists or contained in
  ##     s if it exists
  
  if(missing(selection))
    return("'selection' argument must be supplied")

  ErrMsg <- NULL
  if(length(selection) > 0L && any(is.na(selection)))
    ErrMsg <- c(ErrMsg,
                "'selection' may not be NA")

  if(length(selection) != 1L)
    ErrMsg <- c(ErrMsg,
                "'selection' argument must be a single element vector")

  if(length(selection) > 0L) {
    if(!missing(empty.principal.stratum)) {
      if(!all((selection %in% empty.principal.stratum) | is.na(selection)))
        ErrMsg <- c(ErrMsg,
                    "'selection' value does not appear in specified levels of 's' from 'empty.principal.stratum'")

    } else if(!missing(s)) {
      if(!all((selection %in% s) | is.na(selection)))
        ErrMsg <- c(ErrMsg,
                    "'selection' value does not appear in specified levels of 's' from 'empty.principal.stratum'")
      
    }
  }
  
  return(ErrMsg)
}

.CheckGroupings <- function(groupings, ...) {
  ## groupings must be
  ## 1. not missing
  ## 2. cannot contain NAs
  
  if(missing(groupings))
    return("'groupings' argument must be supplied")

  ErrMsg <- NULL
  if(is.null(groupings) || length(groupings) != 2L)
    ErrMsg <- c(ErrMsg,
                "'groupings' argument must be a two element vector")

  if(length(groupings) && any(is.na(groupings)))
    ErrMsg <- c(ErrMsg,
                "'groupings' may not contain a NA")
  
  return(ErrMsg)
}

.CheckTrigger <- function(trigger, d) {
  ## trigger must be
  ## 1. not missing
  ## 2. cannot contain NAs

  if(missing(trigger))
    return("'trigger' argument must be supplied")

  ErrMsg <- NULL
  if(is.null(trigger) || length(trigger) < 1L)
    ErrMsg <- c(ErrMsg,
                "'trigger' argument must be a single element vector")

  if(length(trigger) && any(is.na(trigger)))
    ErrMsg <- c(ErrMsg,
                "'trigger' may be NA")

  if(length(trigger) > 0L && !missing(d) &&
     !all((trigger %in% d) | is.na(trigger)))
        ErrMsg <- c(ErrMsg,
                    "'trigger' value does not appear in given values of 'd'")
  return(ErrMsg)
}

.CheckLength <- function(z, s, d, y) {
  ## for vectors z, s, d, and y
  ## 1. must be same length

  vectorMissing <- c(missing(z), missing(s), missing(d), missing(y))
  
  if(!any(vectorMissing) &&
     (length(z) != length(s) || length(z) != length(d) ||
      length(z) != length(y)))
    return(sprintf("Vectors %s are of unequal lengths", paste(c('z','s','d','y')[!vectorMissing], collapse=',')))

  return(NULL)
}

.CheckZ <- function(z, groupings, na.rm, ...) {
  ## z must be
  ## 1. not missing
  ## 2. cannot contain NAs unless na.rm == TRUE
  ## 3. all z values must match the values of groupings unless groupings is
  ##      missing

  if(missing(z))
    return("'z' argument must be supplied")

  ErrMsg <- NULL
  if(any(is.na(z)) && na.rm)
    ErrMsg <- c(ErrMsg,
                "'z' cannot contain any NA values")

  if(!missing(groupings) && !all(z %in% groupings | is.na(z)))
    ErrMsg <- c(ErrMsg,
                "All values of 'z' must match one of the two values in 'groupings'")
  return(ErrMsg)
}

.CheckS <- function(s, empty.principal.stratum, na.rm, ...) {
  ## s must be
  ## 1. not missing
  ## 2. cannot contain NAs unless na.rm == TRUE
  ## 3. all s values must match the values of empty.principal.stratum if exists

  if(missing(s))
    return("'s' argument must be supplied")

  ErrMsg <- NULL
  if(any(is.na(s)) && !na.rm)
    ErrMsg <- c(ErrMsg,
                "argument 's' cannot contain any NA values")

  if(!missing(empty.principal.stratum) &&
     !all(s %in% empty.principal.stratum | is.na(s)))
    ErrMsg <- c(ErrMsg,
                "All values of 's' must match one of the two values in 'empty.principal.stratum'")

  return(ErrMsg)
}



.CheckY <- function(y, s, selection, ...) {
  ## y must be
  ## 1. not missing
  ## 2. cannot be NA is s is selected.

  if(missing(y))
    return("'y' argument must be supplied")

  ErrMsg <- NULL
  if(!missing(selection) && !missing(s) && length(s) == length(y) &&
     any(selection %in% s) &&
     any(s %in% selection & !is.na(s) & is.na(y)))
    ErrMsg <- c(ErrMsg,
                sprintf("argument 'y' cannont contain a NA value if the corresponding 's' is %s",
                        paste(selection, collapse=",")))
  return(ErrMsg)
}

.CheckD <- function(d, s, selection) {
  ## d must be
  ## 1. not missing
  ## 2. cannot be NA is s is selected.

  if(missing(d))
    return("'d' argument must be supplied")

  ErrMsg <- NULL
  if(!missing(selection) && !missing(s) && length(s) == length(d) &&
     any(selection %in% s) &&
     any(s %in% selection & !is.na(s) & is.na(d)))
    ErrMsg <- c(ErrMsg,
                sprintf("argument 'd' cannot contain a NA value if the corresponding 's' is %s", paste(selection, collapse=",")))

  return(ErrMsg)
}

.ComputePiPsiPhi <- function(call = match.call(definition=sys.function(sys.parent()),
                               call=sys.call(sys.parent())),
                             envir=parent.frame(n=2),
                             p0, p1) {
  ## Assume only one of Pi, psi, or phi is present.  Calculate other two values.
  .ComputeFunc <- function(psi, Pi, phi, ..., p0, p1){
    if(!missing(psi) && !is.null(psi)) {
      Pi <-
        ifelse(abs(psi) < sqrt(.Machine$double.eps), p0*p1,
               -(sqrt((p1^2-2*p0*p1+p0^2)*exp(2*psi)+p1^2
                      +exp(psi)
                      *(-2*p1^2+2*p1-2*p0^2+2*p0)
                      +(2*p0-2)*p1+p0^2-2*p0+1)
                 +p1+exp(psi)*(-p1-p0)+p0-1)
               /(2*exp(psi)-2))

      phi <- Pi/p1
    } else if(!missing(phi) && !is.null(phi)) {
      Pi <- p1*phi
      psi <- log((p1 * phi^2 + (1 - p0 - p1)*phi)/
                 (p1 * phi^2 - (p1 + p0)* phi + p0))
    } else {
      psi <- log(Pi * (1 - p1 - p0 + Pi)/(p1 - Pi)/(p0 - Pi))
      phi <- Pi/p1
    }

    return(list(psi = psi, Pi=Pi, phi=phi))
  }

  call[[1]] <- .ComputeFunc
  eval(call, envir=data.frame(p0=p0, p1=p1), enclos=envir)
}

.RunCheck <- function(checks,
                      parentCall=match.call(definition=sys.function(sys.parent()),
                        call=sys.call(sys.parent())),
                      envir=parent.frame(n=2)) {

  ErrMsg <- NULL
  for(check in checks) {
    if(is.character(check)) {
      check <- as.symbol(check)
    }

    parentCall[[1]] <- check
    ErrMsg <- c(ErrMsg, eval(parentCall, envir=envir))
  }

  return(ErrMsg)
}

.RunCompute <- function(check, ...,
                        parentCall=match.call(definition=sys.function(sys.parent()),
                          call=sys.call(sys.parent())),
                        envir=parent.frame(n=2)) {
    pNames <- names(parentCall)

    if(is.null(pNames))
      pNames <- ""
    
    thisCall <- match.call(expand.dots=FALSE)

    dotNames <- names(thisCall[['...']])
    dotNames <- dotNames[dotNames != ""]
    
    ## Ensure no name overlap
    if(any(dotNames  %in% names(parentCall))) {
      parentCall[names(thisCall[['...']]) %in% names(parentCall)] <- NULL
    }

    parentCall[seq.int(from=length(parentCall)+1L,
                       length.out=length(dotNames))] <- match.call()[dotNames]
    names(parentCall) <- c(pNames, dotNames)

    parentCall[[1]] <- as.name(check)
    eval(parentCall, envir=parent.frame())
}
