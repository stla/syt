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
#' @seealso \code{\link{KostkaNumber}}.
#'
#' @examples
#' count_ssytx(c(4, 3, 3, 2), 5)
count_ssytx <- function(lambda, n) {
  lambda <- checkPartition(lambda)
  n <- as.integer(n)
  stopifnot(n >= 1L)
  if(n < length(lambda)) {
    return(0L)
  }
  if(n == 1L) {
    return(1L)
  }
  lambda <- c(lambda, rep(0, n-length(lambda)))
  out <- 1
  for(j in 2L:n) {
    for(i in 1L:(j-1L)) {
      out <- out * (1 + (lambda[i] - lambda[j]) / (j - i)) 
    }
  }
  as.integer(ceiling(out))
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
#' @seealso \code{\link{ssytx_withGivenShapeAndWeight}}.
#'
#' @examples
#' ssytx <- all_ssytx(c(2, 1), 3)
#' lapply(ssytx, prettyTableau)
all_ssytx <- function(lambda, n) {
  lambda <- checkPartition(lambda)
  stopifnot(isPositiveInteger(n))
  row <- function(n, len, prev, xxs) {
    if(len == 0L) {
      list(integer(0L))
    } else {
      x <- xxs[[1L]]
      xs <- xxs[-1L]
      do.call(c, lapply(.rg(max(x, prev), n), function(a) {
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

#' @title Semistandard Young tableaux with given shape and weight
#' @description Enumeration of all semistandard Young tableaux with a given 
#'   shape and a given weight. The \emph{weight} of a tableau is the 
#'   vector whose \eqn{i}-th element is the number of occurrences of \eqn{i} 
#'   in this tableau.
#' 
#' @param lambda integer partition, the shape
#' @param weight integer vector, the weight
#'
#' @return List of all semistandard Young tableaux with shape \code{lambda} 
#'   and weight \code{weight}.
#' @export
#' @seealso \code{\link{all_ssytx}}.
#'
#' @examples
#' ssytx <- ssytx_withGivenShapeAndWeight(c(4, 1), c(0, 2, 1, 1, 1))
#' lapply(ssytx, prettyTableau)
ssytx_withGivenShapeAndWeight <- function(lambda, weight) {
  GTpatterns <- GelfandTsetlinPatterns(lambda, weight)
  lapply(GTpatterns, .GTpatternToTableau)
}