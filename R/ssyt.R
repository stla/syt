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
  out
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
#' all_ssytx(c(2, 1), 3)
all_ssytx <- function(lambda, n) {
  lambda <- as.integer(checkPartition(lambda))
  n <- as.integer(n)
  stopifnot(n >= 1L)
  l <- length(lambda)
  ll <- lambda[l]
  tmp <- list(list())
  res <- list()
  while(length(tmp) > 0L) {
    T <- tmp[[1L]]
    tmp <- tmp[-1L]
    k <- length(T)
    if(k == l && length(T[[k]]) >= ll) {
      res <- c(res, list(T))
    } else {
      if(k == 0L || length(T[[k]]) %in% c(0L, lambda[k])) {
        if(k == 0L) {
          start <- 1L
        } else {
          start <- T[[k]][1L] + 1L
        }
        for(u in rg(start, n)) {
          U <- c(T, list(u))
          tmp <- c(tmp, list(U))
        }
      } else {
        Tk <- T[[k]]
        lk <- length(Tk)
        if(k == 1L) {
          start <- Tk[lk]
        } else {
          start <- max(Tk[lk], T[[k-1L]][lk] + 1L)
        }
        for(u in rg(start, n)) {
          U <- T
          U[[k]] <- c(Tk, u)
          tmp <- c(tmp, list(U))
        }
      }
    }
  }
  res
} 
