.CheckEmptyPrincipleStratum <-function(call=match.call(definition=sys.function(sys.parent()),
                                         call=sys.call(sys.parent()))){
  .CheckFunc <- function(empty.principle.stratum, ...) {
    ## empty.principle.stratum must be
    ## 1. not missing
    ## 2. not NULL
    ## 3. a 2 element vector
    ## 4. elements must be diffrent values
    ## 5. cannot be NA

    if(missing(empty.principle.stratum))
      return("'empty.principle.stratum' argument must be supplied")

    ErrMsg <- NULL
    if(is.null(empty.principle.stratum) ||
       length(empty.principle.stratum) != 2L)
      ErrMsg <- c(ErrMsg,
                  "'empty.principle.stratum' argument must be a two element vector")
    else if(empty.principle.stratum[1] == empty.principle.stratum[2])
      ErrMsg <- c(ErrMsg,
                  "'empty.principle.stratum' argument cannot have the same values for both elements")
    

    if(length(empty.principle.stratum) &&
       any(is.na(empty.principle.stratum)))
      ErrMsg <- c(ErrMsg,
                  "'empty.principle.stratum' may not contain a NA")

    return(ErrMsg)
  }

  call[[1]] <- as.name(".CheckFunc")
  eval(call)
}

.CheckSelection <- function(call = match.call(definition=sys.function(sys.parent()),
                      call=sys.call(sys.parent()))) {
  
  .CheckFunc <- function(selection, s, empty.principle.stratum, ...) {
    ## selection must be
    ## 1. not missing
    ## 2. a single element vector
    ## 3. not NA
    ## 4. contained in empty.principle.stratum if it exists or contained in
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
      if(!missing(empty.principle.stratum)) {
        if(!all((selection %in% empty.principle.stratum) | is.na(selection)))
      ErrMsg <- c(ErrMsg,
                  "'selection' value does not appear in specified levels of 's' from 'empty.principle.stratum'")

      } else if(!missing(s)) {
        if(!all((selection %in% s) | is.na(selection)))
        ErrMsg <- c(ErrMsg,
                    "'selection' value does not appear in specified levels of 's' from 'empty.principle.stratum'")
     
      }
    }
    
    return(ErrMsg)
  }

  call[[1]] <- as.name(".CheckFunc")
  eval(call)
}

.CheckGroupings <- function(call = match.call(definition=sys.function(sys.parent()),
                      call=sys.call(sys.parent()))) {
  .CheckFunc <- function(groupings, ...) {
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
    
    ErrMsg
  }
  call[[1]] <- as.name(".CheckFunc")
  eval(call)
}

.CheckZ <- function(call = match.call(definition=sys.function(sys.parent()),
                      call=sys.call(sys.parent()))) {
  ## z must be
  ## 1. not missing
  ## 2. cannot contain NAs unless na.rm == TRUE
  ## 3. all z values must match the values of groupings unless groupings is missing
  
  .CheckFunc <- function(z, groupings, na.rm, ...) {
    if(missing(z))
      return("'z' argument must be supplied")

    ErrMsg <- NULL
    if(any(is.na(z)) && na.rm)
      ErrMsg <- c(ErrMsg,
                  "'z' cannot contain any NA values")

    if(!missing(groupings) && !all(z %in% groupings | is.na(z)))
      ErrMsg <- c(ErrMsg,
                  "All values of 'z' must match one of the two values in 'groupings'")

  }

  call[[1]] <- as.name(".CheckFunc")
  eval(call)
}
  
.CheckS <- function(call = match.call(definition=sys.function(sys.parent()),
                      call=sys.call(sys.parent()))) {
  ## s must be
  ## 1. not missing
  ## 2. cannot contain NAs unless na.rm == TRUE
  ## 3. all s values must match the values of empty.principle.stratum if exists
  .CheckFunc <- function(s, empty.principle.stratum, ...) {
    if(missing(s))
      return("'s' argument must be supplied")

    ErrMsg <- NULL
    if(any(is.na(s)) && !na.rm)
      ErrMsg <- c(ErrMsg,
                  "argument 's' cannot contain any NA values")

    if(!missing(empty.principle.stratum) &&
       !all(s %in% empty.principle.stratum | is.na(s)))
      ErrMsg <- c(ErrMsg,
                  "All values of 's' must match one of the two values in 'empty.principle.stratum'")

    return(ErrMsg)
  }

  call[[1]] <- as.name(".CheckFunc")
  eval(call)
}  

.CheckY <- function(call = match.call(definition=sys.function(sys.parent()),
                      call=sys.call(sys.parent()))) {
  ## s must be
  ## 1. not missing
  ## 2. cannot be NA is s is selected.
  .CheckFunc <- function(y, s, selection, ...) {
    if(missing(y))
      return("'y' argument must be supplied")

    ErrMsg <- NULL
    if(!missing(selection) && !missing(s) && length(s) == length(y) &&
       any(selection %in% s) &&
       any(s %in% selection & !is.na(s) & is.na(y)))
      ErrMsg <- c(ErrMsg,
                  sprintf("argument 'y' cannont contain a NA value if the corresponding 's' is %s",
                          paste(selection, collapse=",")))
  }

  call[[1]] <- as.name(".CheckFunc")
  eval(call)
}

.ComputePiPsiPhi <- function(call = match.call(definition=sys.function(sys.parent()),
                               call=sys.call(sys.parent())),
                             p0, p1) {
  ## Assume only one of Pi, psi, or phi is present.  Calculate other two values.
  .ComputeFunc <- function(psi, Pi, phi, ...) {
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

  call[[1]] <- as.name(".ComputeFunc")
  eval(call)
}
  
.RunCheck <- function(check, ...) {
    parentCall <- match.call(definition=sys.function(sys.parent()),
                             call=sys.call(sys.parent()))

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
