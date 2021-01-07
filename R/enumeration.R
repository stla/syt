ytbl <- function(N, lambda, a, more){

  #N <- it <- sum(lambda)
  it <- N

  if(more){
    lambda <- integer(N); lambda[1L] <- 1L
    isave <- 0L

    for(i in 2L:N){
      lambda[a[i]] <- lambda[a[i]] + 1L
      if(a[i] < a[i-1L]){
        isave <- i
        break
      }
    }

    if(isave == 0L){
      return(list(a=a, more=FALSE))
    }

    it <- lambda[a[isave]+1L]

    for(i in N:1L){
      if(lambda[i] == it){
        a[isave] <- i
        lambda[i] <- lambda[i] - 1L
        it <- isave - 1L
        break
      }
    }
  }  # end of if(more)

  lambda0 <- integer(N); lambda0[seq_along(lambda)] <- as.integer(lambda)
  lambda <- lambda0
  k <- ir <- 1L

  while(TRUE){
    if(N < ir){
      break
    }
    if(lambda[ir] != 0L){
      a[k] <- ir
      lambda[ir] <- lambda[ir] - 1L
      k <- k+1L
      ir <- ir+1L
      next
    }
    if(it < k){
      break
    }
    ir <- 1L
  }

  if(N == 1L){
    return(list(a=a, more=FALSE))
  }

  for(j in seq_len(N-1L)){
    if(a[j+1L] < a[j]){
      return(list(a=a, more=TRUE))
    }
  }

  return(list(a=a, more=FALSE))
}

#' Enumeration of standard Young tableaux
#' @description Generates all standard Young tableaux of a given shape.
#' @param lambda shape, an integer partition
#'
#' @return A list of standard Young tableaux.
#' @export
#'
#' @examples
#' all_sytx(c(5,2))
all_sytx <- function(lambda){
  lambda <- checkPartition(lambda)
  N <- sum(lambda)
  nextA <- ytbl(N, lambda, integer(N), FALSE)
  As <- list(nextA$a)
  while(nextA$more){
    nextA <- ytbl(N, lambda, nextA$a, TRUE)
    As <- c(As, list(nextA$a))
  }
  lapply(As, .ballot2syt)
}

.islastsyt <- function(syt){
  a <- .syt2ballot(syt)
  N <- length(a)
  for(j in seq_len(N-1L)){
    if(a[j+1L] < a[j]){
      return(FALSE)
    }
  }
  TRUE
}

.firstsyt <- function(lambda, a){
  lambda <- removezeros(lambda)
  N <- it <- sum(lambda)

  if(N==1L){
    return(list(1L))
  }

  lambda0 <- integer(N)
  lambda0[seq_along(lambda)] <- as.integer(lambda)
  lambda <- lambda0
  k <- ir <- 1L

  while(TRUE){
    if(N < ir){
      break
    }
    if(lambda[ir] != 0L){
      a[k] <- ir
      lambda[ir] <- lambda[ir] - 1L
      k <- k+1L
      ir <- ir+1L
      next
    }
    if(it < k){
      break
    }
    ir <- 1L
  }

  .ballot2syt(a)
}

.nextsyt <- function(syt){
  if(.islastsyt(syt)){
    return(NULL)
  }
  a <- .syt2ballot(syt)
  N <- sum(lengths(syt))
  lambda <- integer(N); lambda[1L] <- 1L
  isave <- 0L;

  for(i in 2L:N){
    lambda[a[i]] <- lambda[a[i]] + 1L
    if(a[i] < a[i-1L]){
      isave <- i
      break
    }
  }

  if(isave == 0L){
    return(.ballot2syt(a))
  }

  it <- lambda[a[isave]+1L]

  for(i in N:1L){
    if(lambda[i] == it){
      a[isave] <- i
      lambda[i] <- lambda[i] - 1L
      it <- isave - 1L
      break
    }
  }

  .firstsyt(lambda, a)
}


#' First tableau of a given shape
#' @description Returns the "first" standard Young tableau of a given shape.
#'
#' @param lambda the shape, an integer partition
#'
#' @return A standard Young tableau.
#' @export
#'
#' @examples
#' firstsyt(c(4,2,1))
firstsyt <- function(lambda){
  lambda <- checkPartition(lambda)
  .firstsyt(lambda, integer(sum(lambda)))
}

#' Next tableau
#' @description Given a standard Young tableau, returns the "next" one having
#' the same shape.
#'
#' @param syt a standard Young tableau
#'
#' @return A standard Young tableau of the same shape as \code{syt}, or
#' \code{NULL} if \code{syt} is the last standard Young tableau of this shape.
#' @export
#'
#' @examples
#' syt <- firstsyt(c(4,2,1))
#' nextsyt(syt)
nextsyt <- function(syt){
  checkSYT(syt)
  .nextsyt(syt)
}
