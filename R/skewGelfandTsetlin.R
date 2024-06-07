# assumes lambda clean and mu has trailing zeros in order that length(lambda)=length(mu)
sandwichedPartitions <- function(weight, mu, lambda) {
  # assumes length(a_as) == length(b_bs)
  recursiveFun <- function(d, h0, a_as, b_bs) {
    if(d < 0L || d < sum(a_as) || d > sum(b_bs)) {
      list()
    } else if(d == 0L) {
      list(rep(0L, length(a_as)))
    } else {
      a <- a_as[1L]
      as <- a_as[-1L]
      b <- b_bs[1L]
      bs <- b_bs[-1L]
      hrange <- .rg(max(0L, a), min(h0, b))
      do.call(c, lapply(hrange, function(h) {
        lapply(recursiveFun(d-h, h, as, bs), function(hs) {
          c(h, hs)
        })
      }))
    }
  }
  recursiveFun(weight, lambda[1L], mu, lambda)
}


#' @title Skew Gelfand-Tsetlin patterns
#' @description Enumeration of Gelfand-Tsetlin patterns defined by a 
#'   given skew partition and a given weight.
#' 
#' @param lambda,mu integer partitions defining the skew partition: 
#'   \code{lambda} is the outer partition and \code{mu} is the inner partition 
#'   (so \code{mu} must be a subpartition of \code{lambda}); \code{lambda}
#'   will be the last row of the generated Gelfand-Tsetlin patterns and 
#'   \code{mu} will be their first row
#' @param weight integer vector; this vector will be the 
#'   differences of the row sums of the generated Gelfand-Tsetlin patterns; 
#'   consequently, there will be no generated Gelfand-Tsetlin pattern unless 
#'   the sum of \code{weight} equals the difference between the sum of 
#'   \code{lambda} and the sum of \code{mu}
#'
#' @return A list of matrices with non-negative integer entries. The number 
#'   of columns of these matrices is the length of \code{lambda} and the 
#'   number of rows of these matrices is one plus the length of \code{weight}.
#' @export
#' @importFrom igraph graph_from_edgelist V all_simple_paths
#'
#' @examples
#' skewGelfandTsetlinPatterns(c(3, 1, 1), c(2), c(1, 1, 1))
skewGelfandTsetlinPatterns <- function(lambda, mu, weight) {
  stopifnot(isPartition(lambda), isPartition(mu))
  stopifnot(isIntegerVector(weight))
  lambda <- as.integer(removezeros(lambda))
  mu <- as.integer(removezeros(mu))
  ellLambda <- length(lambda)
  ellMu <- length(mu)
  if(ellLambda < ellMu) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  wMu <- sum(mu)
  mu <- c(mu, rep(0L, ellLambda - ellMu))
  if(any(lambda < mu)) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  if(any(weight < 0L)) {
    return(list())
  }
  wLambda <- sum(lambda)
  wWeight <- sum(weight)
  if(wWeight != wLambda - wMu) {
    return(list())
  }
  if(wWeight == 0L) {
    return(
      list(do.call(
        rbind,
        replicate(length(weight) + 1L, lambda, simplify = FALSE)
      ))
    )
  }
  #
  recursiveFun <- function(kappa, w) {
    d <- sum(kappa) - w[length(w)]
    if(d == wMu) {
      if(all(kappa >= mu) &&
         all(head(mu, -1L) >= tail(kappa, -1L))
      ) {
        return(list(rbind(mu, kappa, deparse.level = 0L)))
      } else {
        return(list())
      }
    }
    partitions <- sandwichedPartitions(d, c(tail(kappa, -1L), 0L), kappa)
    hw <- head(w, -1L)
    do.call(
      c,
      lapply(partitions, function(nu) {
        lapply(recursiveFun(nu, hw), function(M) {
          rbind(M, kappa, deparse.level = 0L)
        })
      })  
    )
  }
  patterns <- recursiveFun(lambda, weight[weight != 0L])
  if(any(weight == 0L)) {
    indices <- cumsum(pmin(1L, c(1L, weight)))  
    patterns <- lapply(patterns, function(pattern) {
      pattern[indices, , drop = FALSE]
    })
  }
  patterns
}

# convert a skew Gelfand-Tsetlin pattern to a semistandard skew tableau
.skewGTpatternToTableau <- function(pattern) {
  if(ncol(pattern) == 0L) {
    return(list())
  }
  mu <- pattern[1L, ]
  skewTableau <- lapply(mu, function(i) {
    rep(NA_integer_, i)
  })
  partitions <- apply(pattern, 1L, removeTrailingZeros, simplify = FALSE)
  for(i in 2L:nrow(pattern)) {
    skewTableau <- 
      .growTableau(i-1L, skewTableau, partitions[[i]], partitions[[i-1L]])
  }
  skewTableau
}