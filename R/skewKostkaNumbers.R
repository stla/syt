insertWith <- function(f, mp, key, value) {
  if(key %in% names(mp)) {
    mp[[key]] <- f(value, mp[[key]])
  } else {
    mp[[key]] <- value
  }
  mp
}

zip3 <- function(v1, v2, v3) {
  lapply(1L:length(v1), function(i) {
    c(v1[i], v2[i], v3[i])
  })
}

fillings <- function(n, diagram) {
  if(nrow(diagram) == 0L) {
    list(list(integer(0L), integer(0L)))
  } else {
    xy <- diagram[1L, ]
    x <- xy[1L]
    y <- xy[2L]
    rest <- diagram[-1L, , drop = FALSE]
    diagram <- apply(diagram, 1L, partitionAsString)
    upper <-
      n + 1L - match(partitionAsString(c(x, y + 1L)), diagram, nomatch = n + 1L)
    lower <-
      n + 1L - match(partitionAsString(c(x - 1L, y)), diagram, nomatch = n + 1L)
    L <- lapply(fillings(n - 1L, rest), function(filling) {
      nextLetter(lower, upper, filling)
    })
    do.call(c, L)
  }
}

nextLetter <- function(lower, upper, filling) {
  nu <- filling[[1L]]
  lpart <- filling[[2L]]
  shape <- c(nu, 0L)
  lb <- if(lower > 0L) lpart[lower] else 0L
  ub <- if(upper > 0L) min(length(shape), lpart[upper]) else length(shape)
  f <- function(j) {
    if(j == 1L || shape[j-1L] > shape[j]) j else 0L
  }
  v <- vapply((lb+1L):ub, f, integer(1L))
  nlist <- v[v > 0L]
  lapply(nlist, function(i) {
    list(incr(i, shape), c(lpart, i))
  })
}

incr <- function(i, xxs) {
  if(length(xxs) == 0L) {
    integer(0L)
  } else {
    if(i == 0L) {
      finish(xxs)
    } else if(i == 1L) {
      c(xxs[1L] + 1L, finish(xxs[-1L]))
    } else {
      c(xxs[1L], incr(i - 1L, xxs[-1L]))
    }
  }
}

finish <- function(xxs) {
  if(length(xxs) == 0L) {
    integer(0L)
  } else {
    x <- xxs[1L]
    if(x > 0L) {
      c(x, finish(xxs[-1L]))
    } else {
      integer(0L)
    }
  }
}

#' @title Littlewood-Richardson rule for skew Schur polynomial
#' @description Expression of a skew Schur polynomial as a linear
#'   combination of Schur polynomials.
#'
#' @param lambda,mu integer partitions defining the skew partition:
#'   \code{lambda} is the outer partition and \code{mu} is the inner partition
#'   (so \code{mu} must be a subpartition of \code{lambda})
#' @param output the type of the output, \code{"dataframe"} or \code{"list"}
#'
#' @return This computes the expression of the skew Schur polynomial
#'   associated to the skew partition defined by \code{lambda} and \code{mu}
#'   as a linear combination of Schur polynomials. If \code{output="dataframe"},
#'   the output is a dataframe with two columns: the column \code{coeff} gives
#'   the coefficients of this linear combination, and the column \code{nu}
#'   gives the partitions defining the Schur polynomials of this linear
#'   combination as character strings, e.g. the partition \code{c(4, 3, 1)} is
#'   given by \code{"4, 3, 1"}. If \code{output="list"}, the output is a list
#'   with two fields: the field \code{coeff} is the vector made of the
#'   coefficients of the linear combination, and the field \code{nu} is the
#'   list of partitions defining the Schur polynomials of the linear combination
#'   given as integer vectors.
#' @noRd
LRskew <- function(lambda, mu, output = "dataframe") {
  output <- match.arg(output, c("list", "dataframe"))
  # lambda <- as.integer(removeTrailingZeros(lambda))
  # mu <- as.integer(removeTrailingZeros(mu))
  # ellLambda <- length(lambda)
  # ellMu <- length(mu)
  # if(ellLambda < ellMu) {
  #   stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  # }
  # mu <- c(mu, rep(0L, ellLambda - ellMu))
  # if(any(lambda < mu)) {
  #   stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  # }
  ellLambda <- length(lambda)
  n <- sum(lambda) - sum(mu)
  if(n == 0L) {
    if(output == "dataframe") {
      return(data.frame("coeff" = 1L, "nu" = "[]"))
    } else {
      return(list("coeff" = 1L, "nu" = list(integer(0L))))
    }
  }
  mu <- c(mu, rep(0L, ellLambda - length(mu)))
  f <- function(old, nu) {
    insertWith(`+`, old, partitionAsString(nu), 1L)
  }
  Liab <- rev(zip3(seq_len(ellLambda), lambda, mu))
  diagram <- do.call(rbind, do.call(c, lapply(Liab, function(iab) {
    i <- iab[1L]
    a <- iab[2L]
    b <- iab[3L]
    jvec <- if(b < a) (b + 1L):a else integer(0L)
    lapply(jvec, function(j) {
      c(i, j)
    })
  })))
  Lnu <- lapply(fillings(n, diagram), `[[`, 1L)
  v <- Reduce(f, Lnu, init = integer(0L))
  if(output == "dataframe") {
    data.frame("coeff" = v, "nu" = names(v))
  } else {
    partitions <- lapply(names(v), fromPartitionAsString)
    list("coeff" = unname(v), "nu" = partitions)
  }
}

lastSubpartition <- function(w, lambda) {
  if(w == 0L || length(lambda) == 0L) {
    integer(0L)
  } else {
    k <- lambda[1L]
    if(w <= k) {
      w
    } else {
      c(k, lastSubpartition(w - k, tail(lambda, -1L)))
    }
  }
}

#' @title Skew Kostka numbers
#' @description Skew Kostka numbers associated to a given skew partition.
#' 
#' @param lambda,mu integer partitions defining the skew partition: 
#'   \code{lambda} is the outer partition and \code{mu} is the inner partition 
#'   (so \code{mu} must be a subpartition of \code{lambda})
#' @param output the format of the output, either \code{"vector"} or 
#'   \code{"list"}
#'
#' @return If \code{output="vector"}, the function returns a named vector. 
#'   This vector is made of the positive skew Kostka numbers 
#'   \eqn{K_{\lambda/\mu,\nu}} and its names encode the partitions \eqn{\nu}.
#'   If \code{ouput="list"}, the function returns a list. Each element of this 
#'   list is a named list with two elements: an integer partition \eqn{\nu} 
#'   in the field named \code{"nu"}, and the corresponding skew Kostka number 
#'   \eqn{K_{\lambda/\mu,\nu}} in the field named \code{"value"}. Only the 
#'   non-null skew Kostka numbers are provided by this list.
#' @export
#' @importFrom partitions parts
#' @seealso \code{\link{KostkaNumber}}, \code{\link{KostkaNumbersWithGivenMu}}.
#'
#' @examples
#' skewKostkaNumbers(c(4,2,2), c(2,2))
skewKostkaNumbers <- function(lambda, mu, output = "vector") {
  stopifnot(isPartition(lambda), isPartition(mu))
  output <- match.arg(output, c("vector", "list"))
  lambda <- as.integer(removeTrailingZeros(lambda))
  mu <- as.integer(removeTrailingZeros(mu))
  ellLambda <- length(lambda)
  ellMu <- length(mu)
  if(ellLambda < ellMu || any(head(lambda, ellMu) < mu)) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  kappa <- lastSubpartition(sum(lambda)-sum(mu), lambda)
  nus <- .dominatedPartitions(kappa)
  lr <- LRskew(lambda, mu, output = "list")
  pis <- lr[["nu"]]
  coeffs <- lr[["coeff"]]
  listOfIndexVectors <- lapply(nus, function(nu) {
    which(vapply(pis, function(pi) {
      .isDominatedBy(nu, pi)
    }, logical(1L), USE.NAMES = FALSE))
  })
  indices <- which(lengths(listOfIndexVectors) != 0L)
  pisAsStrings <- 
    vapply(pis, partitionAsString, character(1L), USE.NAMES = FALSE)
  kNumbers <- vapply(indices, function(j) {
    nu <- nus[[j]]
    kNumbers_nu <- .KostkaNumbersWithGivenMu(nu)
    i_ <- listOfIndexVectors[[j]]
    sum(coeffs[i_] * vapply(pisAsStrings[i_], function(piAsString) {
      kNumbers_nu[[piAsString]][["value"]]
    }, integer(1L), USE.NAMES = FALSE))
  }, integer(1L), USE.NAMES = FALSE)
  nus <- nus[indices]
  names(kNumbers) <-
    vapply(nus, partitionAsString, character(1L), USE.NAMES = FALSE)
  if(output == "vector") {
    kNumbers
  } else {
    mapply(
      function(value, nu) {
        list("nu" = nu, "value" = value)
      },
      kNumbers, nus,
      USE.NAMES = TRUE, SIMPLIFY = FALSE
    )
  }
}
