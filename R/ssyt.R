#' Number of semistandard Young tableaux
#' @description Number of semistandard Young tableaux of a given shape and 
#'   filled with integers between \code{1} and a given \code{n}. 
#' 
#' @param lambda an integer partition, the shape
#' @param n an integer, the maximum value of the entries (the minimum value 
#'   is \code{1})
#'
#' @return The number of semistandard Young tableaux with shape \code{lambda} 
#'   and filled with integers between \code{1} and \code{n}.
#' @export
#'
#' @examples
#' count_ssytx(c(4, 3, 3, 2), 5)
count_ssytx <- function(lambda, n) {
  lambda <- checkPartition(lambda)
  n <- as.integer(n)
  stopifnot(n >= 1L)
  if(n < length(lambda)) {
    return(0)
  }
  if(n == 1L) {
    return(1)
  }
  lambda <- c(lambda, rep(0, n-length(lambda)))
  out <- 1
  for(j in 2L:n) {
    for(i in 1L:(j-1L)) {
      out <- out * (1 + (lambda[i] - lambda[j]) / (j - i)) 
    }
  }
  as.integer(out)
}

rg <- function(start, end) {
  seq_len(end - start + 1L) + (start - 1L)
}

#' Enumeration of semistandard Young tableaux
#' @description Generates all semistandard Young tableaux of a given shape and 
#'   filled with integers between \code{1} and a given \code{n}. 
#' 
#' @param lambda an integer partition, the shape
#' @param n an integer, the maximum value of the entries (the minimum value 
#'   is \code{1})
#'
#' @return List of all semistandard Young tableaux with shape \code{lambda} 
#'   and filled with integers between \code{1} and \code{n}.
#' @export
#'
#' @examples
#' ssytx <- all_ssytx(c(2, 1), 3)
#' lapply(ssytx, prettyTableau)
all_ssytx <- function(lambda, n) {
  row <- function(n, len, prev, xxs) {
    if(len == 0L) {
      list(integer(0L))
    } else {
      x <- xxs[[1L]]
      xs <- xxs[-1L]
      do.call(c, lapply(rg(max(x, prev), n), function(a) {
        lapply(row(n, len - 1L, a, xs), function(as) {
          c(a, as)
        })
      }))
    }
  }
  worker <- function(prevRow, sss) {
    if(length(sss) == 0L) {
      list(list())
    } else {
      s <- sss[[1L]]
      ss <- sss[-1L]
      do.call(c, lapply(row(n, s, 1L, prevRow), function(r) {
        lapply(worker(r + 1L, ss), function(rs) {
          c(list(r), rs)
        })
      }))
    }
  }
  worker(rep(0L, lambda[1L]), lambda)
}
