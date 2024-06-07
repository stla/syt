diffSequence <- function(x) {
  c(diff(-x), x[length(x)])
}

mkSkewTableau <- function(skewTableau) {
  lapply(skewTableau, function(row) {
    offset <- row[[1L]]
    c(rep(NA_integer_, offset), row[[2L]])
  })
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
  if(any(bs < 0L)) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  ds <- diffSequence(as)
  results <- worker(as, bs, ds, rep(1L, bs[1L]))
  lapply(results, mkSkewTableau)
}

#' @title Skew semistandard tableaux with given shape and weight
#' @description Enumeration of all skew semistandard tableaux with a given 
#'   shape and a given weight. The \emph{weight} of a tableau is the 
#'   vector whose \eqn{i}-th element is the number of occurrences of \eqn{i} 
#'   in this tableau.
#' 
#' @param lambda,mu integer partitions defining the skew partition: 
#'   \code{lambda} is the outer partition and \code{mu} is the inner partition 
#'   (so \code{mu} must be a subpartition of \code{lambda})
#' @param weight integer vector, the weight
#'
#' @return List of all skew semistandard tableaux whose shape is the skew 
#'   partition defined by \code{lambda} and \code{mu} and whose weight is 
#'   \code{weight}.
#' @export
#'
#' @examples
#' ssstx <- skewTableauxWithGivenShapeAndWeight(c(3, 1, 1), c(2), c(1, 1, 1))
#' lapply(ssstx, prettySkewTableau)
skewTableauxWithGivenShapeAndWeight <- function(lambda, mu, weight) {
  skewGTpatterns <- 
    skewGelfandTsetlinPatterns(lambda, mu, removeTrailingZeros(weight))
  lapply(skewGTpatterns, .skewGTpatternToTableau)
}

#' @title Check whether a tableau is a skew tableau
#' @description Check whether a tableau is a skew tableau.
#' 
#' @param tableau a tableau
#'
#' @return A Boolean value.
#' @export
#'
#' @examples
#' tbl <- list(c(NA, NA, 1, 1), c(NA, 1), c(1, 2))
#' isSkewTableau(tbl)
isSkewTableau <- function(tableau) {
  isTableau(tableau) && all(vapply(tableau, function(row) {
    nas <- which(is.na(row))
    length(nas) < length(row) && identical(nas, seq_along(nas))
  }, logical(1L)))
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
  stopifnot(isSkewTableau(skewTableau))
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

#' @title Dual skew tableau
#' @description Returns the dual (skew) tableau of a skew tableau.
#'
#' @param skewTableau a skew tableau
#'
#' @return A skew tableau.
#' @export
#'
#' @examples
#' tbl <- list(c(NA, NA, 1, 1), c(NA, 1), c(1, 2))
#' dtbl <- dualSkewTableau(tbl)
#' prettySkewTableau(dtbl)
dualSkewTableau <- function(skewTableau) {
  stopifnot(isSkewTableau(skewTableau))
  sub <- function(b, sktbl) {
    if(length(sktbl) == 0L) {
      list(b, integer(0L))
    } else {
      athis <- sktbl[[1L]]
      rest <- sktbl[-1L]
      a <- athis[[1L]]
      this <- athis[[2L]]
      if(a > 0L) {
        sub(b + 1L, rest)
      } else {
        L <- c(list(this), lapply(rest, `[[`, 2L))
        L <- Filter(function(x) length(x) > 0L, L)
        ys <- vapply(L, `[`, integer(1L), 1L)
        list(b, ys)
      }
    }
  }
  strip <- function(sktbl) {
    if(length(sktbl) == 0L) {
      list()
    } else {
      axs <- sktbl[[1L]]
      a <- axs[[1L]]
      xs <- axs[[2L]]
      rest <- sktbl[-1L]
      if(a > 0L) {
        c(list(list(a - 1L, xs)), strip(rest))
      } else {
        if(length(xs) == 0L) {
          list()
        } else {
          if(length(xs) == 1L) {
            list()
          } else {
            c(list(list(0L, xs[-1L])), strip(rest))
          }
        }
      }
    }
  }
  go <- function(axs) {
    if(length(axs) == 0L) {
      list()
    } else {
      this <- sub(0L, axs)
      if(this[[1L]] == 0L && length(this[[2L]]) == 0L) {
        list()
      } else {
        c(list(this), go(strip(axs)))
      }
    }
  }
  offsets <- lapply(skewTableau, function(row) {
    sum(is.na(row))
  })
  contents <- lapply(skewTableau, function(row) {
    Filter(Negate(is.na), as.integer(row))
  })
  axs <- lapply(seq_along(offsets), function(i) {
    list(offsets[[i]], contents[[i]])
  })
  mkSkewTableau(go(axs))
}

#' @title Check whether a skew tableau is semistandard
#' @description Check whether a skew tableau is a semistandard skew tableau.
#'
#' @param skewTableau a skew tableau
#'
#' @return A Boolean value.
#' @export
#'
#' @examples
#' tbl <- list(c(NA, NA, 1, 1), c(NA, 1), c(1, 2))
#' isSemistandardSkewTableau(tbl)
isSemistandardSkewTableau <- function(skewTableau) {
  stopifnot(isSkewTableau(skewTableau))
  contents <- Filter(Negate(is.na), unlist(skewTableau))
  checkContents <- arePositiveIntegers(contents)
  checkRows <- vapply(skewTableau, function(row) {
    isWeaklyIncreasing(Filter(Negate(is.na), row))
  }, logical(1L))
  checkColumns <- vapply(dualSkewTableau(skewTableau), function(row) {
    isStrictlyIncreasing(Filter(Negate(is.na), row))
  }, logical(1L))
  checkContents && all(checkRows) && all(checkColumns)
}

#' @title Check whether a skew tableau is standard
#' @description Check whether a skew tableau is a standard skew tableau.
#'
#' @param skewTableau a skew tableau
#'
#' @return A Boolean value.
#' @export
#'
#' @examples
#' tbl <- list(c(NA, NA, 1, 1), c(NA, 1), c(1, 2))
#' isStandardSkewTableau(tbl)
isStandardSkewTableau <- function(skewTableau) {
  stopifnot(isSkewTableau(skewTableau))
  contents <- Filter(Negate(is.na), unlist(skewTableau))
  isSemistandardSkewTableau(skewTableau) && 
    setequal(contents, seq_along(contents))
}
