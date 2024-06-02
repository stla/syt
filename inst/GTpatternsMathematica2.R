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

partitionsFittingRectangleWithZeros <- function(h, w, d) {
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
      lapply(partitionsFittingRectangleWithZeros(i, w-1L, d-i), function(p) {
        c(i, p)
      })
    })
  )
}

findGTpatterns2 <- function(l, mu, w) {
  w1 <- cumsum(rev(w))
  most_w1 <- head(w1, -1L)
  n <- sum(l)
  ell <- length(l)
  m <- max(l)
  restrParts <- lapply(most_w1, function(k) {
    partitionsFittingRectangleWithZeros(m, ell, n - k)
    # rparts(n-k, ell, m)
  })
  partitions <- c(list(list(l)), restrParts, list(list(mu)))
  Tuples <- function(set1, set2) {
    Grid <- as.matrix(expand.grid(seq_along(set1), seq_along(set2)))
    apply(Grid, 1L, function(ij) {
      list(set1[[ij[1L]]], set2[[ij[2L]]])
    }, simplify = FALSE)
  }
  cfs <- do.call(
    c,
    lapply(seq_len(length(partitions)-1L), function(i) {
      Tuples(partitions[[i]], partitions[[i+1L]])
    })
  )
  # Grid <- as.matrix(expand.grid(lapply(lengths(restrParts), seq_len)))
  # # paths <- apply(Grid, 1L, function(combo) {
  # #   lapply(seq_len(combo), function(i) {
  # #     restrParts[[i]][[combo[i]]]
  # #   })
  # # }, simplify = FALSE)
  # # Grid <- cbind(1L, Grid, 1L)
  # pairs <- cbind(
  #   head(seq_len(2L + ncol(Grid)), -1L), 
  #   tail(seq_len(2L + ncol(Grid)), -1L)
  # )
  # cfs <- apply(Grid, 1L, function(combo) {
  #   path <- c(
  #     list(l),
  #     lapply(seq_along(combo), function(i) {
  #       restrParts[[i]][[combo[i]]]
  #     }),
  #     list(mu)
  #   )
  #   apply(pairs, 1L, function(ij) {
  #     path[ij]
  #   }, simplify = FALSE)
  # }, simplify = FALSE)
  # cfs <- unlist(cfs, recursive = FALSE)
  # ip <- lapply(paths, function(path) {
  #   c(list(l), path, list(mu))
  # })
  # tuples <- function(i) {
  #   Grid <- as.matrix(expand.grid(ip[[i]], ip[[i+1L]]))
  #   apply(Grid, 1L, identity, simplify = FALSE)
  # }
  # List <- do.call(
  #   c,
  #   lapply(seq_len(length(ip)-1L), tuples)
  # )
  # List <- lapply(seq_len(length(ip)-1L), tuples)
  # cfs <- unlist(List, recursive = FALSE, use.names = FALSE)
  print("******************")
  print(str(cfs[[10]]))
  print("******************")
  condition <- function(cf) {
    min(cf[[1L]] - cf[[2L]]) >= 0 &&
      min(head(cf[[2]], -1L) - tail(cf[[1]], -1L)) >=0
  }
  trues <- which(vapply(cfs, condition, logical(1L)))
  edgeList <- do.call(rbind, lapply(cfs[trues], function(cf) {
    c(
      paste0(as.character(cf[[1L]]), collapse = "-"),
      paste0(as.character(cf[[2L]]), collapse = "-")
    )
    #cbind(as.character(cf[[1L]]), as.character(cf[[2L]]))
  }))
  print(vapply(cfs[trues], function(cf) as.integer(sum(cf[[1]]-cf[[2]])), integer(1)))
  print(head(edgeList, 20L))
  gr <- graph_from_edgelist(edgeList)
  print(V(gr))
  all_simple_paths(
    gr, 
    from = paste0(as.character(l), collapse = "-"), 
    to = paste0(as.character(mu), collapse = "-"), 
    mode = "out"
  )
}

l <- c(4, 3, 3, 2, 1, 1, 0, 0)
mu <- c(2, 2, 1, 0, 0, 0, 0, 0)
w <- c(3, 3, 2, 1)

findGTpatterns2(l, mu, w)
