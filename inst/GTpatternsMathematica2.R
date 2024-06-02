library(partitions)
library(igraph)

partitionsFittingRectangle <- function(h, w, d) {
  if(d == 0L) {
    return(list(integer(0L)))
  }
  if(h == 0L || w == 0L) {
    if(d == 0L) {
      return(list(integer(0L)))
    } else {
      return(list())
    }
  }
  do.call(
    c,
    lapply(1L:min(d, h), function(i) {
      lapply(partitionsFittingRectangle(i, w-1L, d-i), function(p) {
        c(i, p)
      })
    })
  )
}

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

#' @importFrom igraph graph_from_edgelist V all_simple_paths
# lambda and mu are clean partitions
findGTpatterns2 <- function(lambda, mu, w) {
  ellLambda <- length(lambda)
  ellMu <- length(mu)
  if(ellLambda < ellMu) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  mu <- c(mu, rep(0L, ellLambda - ellMu))
  if(any(lambda < mu)) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  wLambda <- sum(lambda)
  if(sum(w) != wLambda - sum(mu)) {
    return(list())
  }
  rw <- rev(w)
  # in case w contains some zeros - this will be used at the end
  lines <- cumsum(pmin(1L, c(1L, rw)))
  #
  rw <- rw[rw != 0L]
  m <- lambda[1L]
  listsOfPartitions <- lapply(head(cumsum(rw), -1L), function(k) {
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
  edgeList <- do.call(rbind, lapply(edges, function(edge) {
    c(
      toString(edge[[1L]]),
      toString(edge[[2L]])
    )
  }))
  gr <- graph_from_edgelist(edgeList)
  vertices <- t(vapply(
    names(V(gr)), 
    qspray:::fromString, 
    integer(ellLambda),
    USE.NAMES = FALSE)
  )
  paths <- all_simple_paths(
    gr, 
    from = toString(lambda), 
    to = toString(mu), 
    mode = "out"
  )
  lapply(paths, function(path) {
    vertices[path, ][lines, ]
  })
}

lambda <- c(4, 3, 3, 2, 1, 1)
mu <- c(2, 2, 1)
w <- c(3, 3, 2, 1)

gts <- findGTpatterns2(lambda, mu, w)
gts
#str(gts[[1]])
