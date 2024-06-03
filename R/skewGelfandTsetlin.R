.partitionsFittingRectangleWithZeros <- function(h, w, d) {
  if(d == 0L) {
    return(list(rep(0L, w)))
  }
  if(h == 0L || w == 0L) {
    if(d == 0L) {
      return(list(rep(0L, w)))
    } else {
      return(list())
    }
  }
  do.call(
    c,
    lapply(1L:min(d, h), function(i) {
      lapply(.partitionsFittingRectangleWithZeros(i, w-1L, d-i), function(p) {
        c(i, p)
      })
    })
  )
}

.Pairs <- function(set1, set2) {
  Grid <- as.matrix(expand.grid(seq_along(set1), seq_along(set2)))
  apply(Grid, 1L, function(ij) {
    list(set1[[ij[1L]]], set2[[ij[2L]]])
  }, simplify = FALSE)
}

#' @title Gelfand-Tsetlin patterns of a skew partition
#' @description Enumeration of Gelfand-Tsetlin patterns associated to a 
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
  stopifnot(isIntegerVector(weight), length(weight) >= 1L)
  lambda <- as.integer(removezeros(lambda))
  mu <- as.integer(removezeros(mu))
  ellLambda <- length(lambda)
  ellMu <- length(mu)
  if(ellLambda < ellMu) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  mu <- c(mu, rep(0L, ellLambda - ellMu))
  if(any(lambda < mu)) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  if(any(weight < 0L)) {
    return(list())
  }
  wLambda <- sum(lambda)
  if(sum(weight) != wLambda - sum(mu)) {
    return(list())
  }
  if(all(lambda == mu)) {
    return(
      list(rbind(lambda, lambda))
    )
  }
  rweight <- rev(weight)
  # in case weight contains some zeros - this will be used at the end
  lines <- rev(cumsum(pmin(1L, c(1L, rweight))))
  #
  rweight <- rweight[rweight != 0L]
  m <- lambda[1L]
  listsOfPartitions <- lapply(head(cumsum(rweight), -1L), function(k) {
    .partitionsFittingRectangleWithZeros(m, ellLambda, wLambda - k)
  })
  listsOfPartitions <- c(
    list(list(lambda)), 
    listsOfPartitions, 
    list(list(mu))
  )
  potentialEdges <- do.call(
    c,
    lapply(seq_len(length(listsOfPartitions)-1L), function(i) {
      .Pairs(
        listsOfPartitions[[i]], 
        listsOfPartitions[[i+1L]]
      )
    })
  )
  edges <- Filter(
    function(edge) {
      all(edge[[1L]] >= edge[[2L]]) &&
        all(head(edge[[2L]], -1L) >= tail(edge[[1L]], -1L))
    }, 
    potentialEdges
  )
  if(length(edges) == 0L) {
    return(list())
  }
  edgeList <- do.call(rbind, lapply(edges, function(edge) {
    c(
      toString(edge[[1L]]),
      toString(edge[[2L]])
    )
  }))
  gr <- graph_from_edgelist(edgeList)
  vertices <- t(rbind(vapply(
    names(V(gr)), 
    fromString, 
    integer(ellLambda),
    USE.NAMES = FALSE)
  ))
  paths <- all_simple_paths(
    gr, 
    from = toString(lambda), 
    to = toString(mu), 
    mode = "out"
  )
  lapply(paths, function(path) {
    vertices[path, , drop = FALSE][lines, , drop = FALSE]
  })
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