#' Hooks
#' @description Hooks of a given integer partition.
#'
#' @param lambda integer partition
#'
#' @return The hooks of the partition in a list.
#' @seealso \code{\link{hooklengths}}
#' @export
#' @importFrom partitions conjugate
#'
#' @examples
#' hooks(c(4,2))
hooks <- function(lambda){
  lambda <- as.integer(checkPartition(lambda))
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
#' @description Hook lengths of a given integer partition.
#'
#' @param lambda an integer partition
#'
#' @return The hook lengths of the partition, given in a list.
#' @seealso \code{\link{hooks}}
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
#' @param lambda an integer partition, the shape
#'
#' @return An integer, the number of standard Young tableaux of shape \code{lambda}.
#' @seealso \code{\link{all_sytx}}
#' @export
#'
#' @examples
#' count_sytx(c(5,4,1))
#' length(all_sytx(c(5,4,1)))
count_sytx <- function(lambda){
  numterms <- c()
  denterms <- unlist(hooklengths(lambda))
  denterms <- denterms[-which(denterms==1L)]
  for(i in seq_len(sum(lambda))[-1L]){
    if(!is.na(idx <- match(i,denterms))){
      denterms <- denterms[-idx]
    }else{
      numterms <- c(numterms,i)
    }
  }
  prod(numterms)/prod(denterms)
}

