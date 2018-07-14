#' Hooks
#' @description Hooks of a given partition.
#'
#' @param lambda partition
#'
#' @return The hooks of the partition in a list.
#' @export
#' @importFrom partitions conjugate
#'
#' @examples
#' hooks(c(4,2))
hooks <- function(lambda){
  lambda <- as.integer(lambda)
  dlambda <- partitions::conjugate(lambda)
  out <- vector("list", length(lambda))
  for(i in seq_along(lambda)){
    out[[i]] <- vector("list", lambda[i])
    for(j in seq_len(lambda[i])){
      out[[i]][[j]] <- c(NA_integer_, lambda[i]-j+1L)
    }
  }
  for(j in seq_along(dlambda)){
    for(i in seq_len(dlambda[j])){
      out[[i]][[j]][1L] <- dlambda[j]-i+1L
    }
  }
  return(out)
}

#' Hook lengths
#' @description Hook lengths of a given partition.
#'
#' @param lambda a partition
#'
#' @return The hook lengths of the partition, given in a list.
#' @export
#'
#' @examples
#' hooklengths(c(4,2))
hooklengths <- function(lambda){
  h <- hooks(lambda)
  lapply(h, function(x) unlist(lapply(x, function(y) sum(y)-1L)))
}


#' Number of standard Young tableaux
#' @description Number of standard Young tableaux of a given shape.
#'
#' @param lambda a partition, the shape
#'
#' @return The number of standard Young tableaux of shape \code{lambda}.
#' @export
#'
#' @examples
#' count_sytx(c(5,4,1))
#' length(sytx(c(5,4,1)))
count_sytx <- function(lambda){
  N <- sum(lambda)
  factorial(N)/prod(unlist(hooklengths(lambda)))
}

