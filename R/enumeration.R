ytbl <- function(lambda, a, more){

  N <- it <- sum(lambda)

  if(more){
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
    return(list(a=a, more=FALSE));
  }

  for (j in seq_len(N-1L)){
    if(a[j+1L] < a[j]){
      return(list(a=a, more=TRUE))
    }
  }

  return(list(a=a, more=FALSE))
}


#' Enumeration of standard Young tableaux
#' @description Generates all standard Young tableaux of a given shape.
#' @param lambda shape
#'
#' @return A list of standard Young tableaux.
#' @export
#'
#' @examples
#' sytx(c(5,2))
sytx <- function(lambda){
  nextA <- ytbl(lambda, integer(sum(lambda)), FALSE)
  As <- list(nextA$a)
  while(nextA$more){
    nextA <- ytbl(lambda, nextA$a, TRUE)
    As <- c(As, list(nextA$a))
  }
  lapply(As, vector2syt)
}
