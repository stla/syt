diffSequence <- function(x) {
  c(diff(-x), x[length(x)])
}

#' @title Semistandard skew tableaux
#' @description Enumeration of all semistandard skew tableaux with given shape 
#'   and given maximum entry.
#' 
#' @param lambda,mu integer partitions defining the skew partition: 
#'   \code{lambda} is the outer partition and \code{mu} is the inner partition 
#'   (so \code{mu} must be a subpartition of \code{lambda})
#' @param n a positive integer, the maximum entry of the skew tableaux
#'
#' @return The list of all semistandard skew tableaux whose shape is the skew 
#'   partition defined by \code{lambda} and \code{mu} and with maximum entry 
#'   \code{n}.
#' @export
#'
#' @examples
#' ssstx <- all_ssSkewTableaux(c(4, 3, 1), c(2, 2), 2)
#' lapply(ssstx, prettySkewTableau)
all_ssSkewTableaux <- function(lambda, mu, n) {
  stopifnot(isPartition(lambda), isPartition(mu))
  stopifnot(isPositiveInteger(n))
  n <- as.integer(n)
  row <- function(k, m, aas) {
    if(k == 0L) {
      list(integer(0L))
    } else if(length(aas) == 0L) {
      list()
    } else {
      a <- aas[1L]
      as <- aas[-1L]
      maxam <- max(a, m)
      x_ <- if(maxam <= n) maxam:n else integer(0L)
      do.call(c, lapply(x_, function(x) {
        lapply(row(k - 1L, x, as), function(xs) {
          c(x, xs)
        })
      }))
    }
  }
  worker <- function(aas, bbs, dds, lb) {
    if(length(aas) == 0L) {
      list(list())
    } else {
      a <- aas[1L]
      as <- aas[-1L]
      b <- bbs[1L]
      bs <- bbs[-1L]
      d <- dds[1L]
      ds <- dds[-1L]
      do.call(c, lapply(row(b, 1L, lb), function(this) {
        lbprime <- c(rep(1L, d), this + 1L)
        lapply(worker(as, bs, ds, lbprime), function(rest) {
          c(list(list(a, this)), rest)
        })
      }))
    }
  }
  as <- c(as.integer(mu), rep(0L, length(lambda) - length(mu)))
  bs <- as.integer(lambda) - as
  if(any(bs) < 0L) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  ds <- diffSequence(as)
  results <- worker(as, bs, ds, rep(1L, bs[1L]))
  lapply(results, function(skewpart) {
    lapply(skewpart, function(row) {
      offset <- row[[1L]]
      c(rep(NA_integer_, offset), row[[2L]])
    })
  })
}

#' @importFrom Matrix sparseMatrix
#' @noRd
.skewTableau2matrix <- function(tableau) {
  ls <- lengths(tableau)
  nrows <- length(tableau)
  M <- matrix(NA_integer_, nrow = nrows, ncol = ls[1L])
  for(i in seq_len(nrows)) {
    M[i, seq_len(ls[i])] <- tableau[[i]]
  }
  ij <- which(!is.na(M), arr.ind = TRUE)
  sparseMatrix(i = ij[, 1L],
               j = ij[, 2L],
               x = M[ij])
}

#' @title Pretty skew tableau
#' @description Pretty form of a skew tableau.
#'
#' @param skewTableau a skew tableau
#'
#' @return A '\code{noquote}' character matrix.
#' @export
#' @importFrom Matrix formatSparseM
#'
#' @examples
#' tbl <- list(c(NA, NA, 1, 1), c(NA, 1), c(1, 2))
#' prettySkewTableau(tbl)
prettySkewTableau <- function(skewTableau) {
  M <- .skewTableau2matrix(skewTableau)
  ls <- lengths(skewTableau)
  formattedM <- formatSparseM(M)
  n <- ncol(formattedM)
  furtherFormattedM <- t(vapply(seq_along(ls), function(i) {
    row <- formattedM[i, ]
    row[-seq_len(ls[i])] <- ""
    row
  }, character(n)))
  rownames(furtherFormattedM) <- paste0(seq_len(nrow(formattedM)), " ->")
  noquote(furtherFormattedM)
}
