#' Kostka number
#' @description Computes a Kostka number.
#'
#' @param lambda an integer partition
#' @param mu an integer vector whose sum equals the weight (i.e. the sum) of 
#'   \code{lambda}
#'
#' @return The Kostka number corresponding to \code{lambda} and \code{mu}.
#' @export
#' @importFrom utils head tail
#' @seealso \code{\link{KostkaNumbersWithGivenMu}}, 
#'   \code{\link{KostkaNumbersWithGivenLambda}}, 
#'   \code{\link{skewKostkaNumbers}}.
#'
#' @details
#' The Kostka number \eqn{K(\lambda,\mu)} is the number of semistandard 
#'   Young tableaux with shape \eqn{\lambda} and weight \eqn{\mu}. It does 
#'   not depend on the order of the elements of \eqn{\mu} (so one can always 
#'   take an integer partition for \eqn{\mu}). The \emph{weight} is the 
#'   vector whose \eqn{i}-th element is the number of occurrences of \eqn{i} 
#'   in the tableau.
#'   
#' @examples
#' KostkaNumber(c(3,2), c(1,1,1,2))
#' KostkaNumber(c(3,2), c(1,1,2,1))
#' KostkaNumber(c(3,2), c(1,2,1,1))
#' KostkaNumber(c(3,2), c(2,1,1,1))
#' lambda <- c(4, 3, 1)
#' mu <- rep(1, sum(lambda))
#' KostkaNumber(lambda, mu) == count_sytx(lambda) # should be TRUE
KostkaNumber <- function(lambda, mu) {
  lambda <- checkPartition(lambda)
  mu     <- as.integer(mu)
  wmu  <- sum(mu)
  wlam <- sum(lambda)
  if(wlam == 0L) {
    return(as.integer(wmu == 0L))
  }
  if(wmu != wlam) { #|| !jack:::isDominated(mu, lambda)) {
    return(0L)
  }
  nlam <- length(lambda)
  nmu  <- length(mu)
  n <- max(nlam, nmu)
  lambda <- c(lambda, rep(0L, n - nlam))
  mu     <- c(mu, rep(0L, n - nmu))
  revlam <- rev(lambda)
  boundedNonIncrSeqs <- function(h0, aas, bbs) {
    if(length(aas) == 0L || length(bbs) == 0L) {
      list(integer(0L))
    } else {
      a <- aas[1L]
      b <- bbs[1L]
      as <- aas[-1L]
      bs <- bbs[-1L]
      h_ <- .rg(max(0L, a), min(h0, b))
      do.call(c, lapply(h_, function(h) {
        lapply(boundedNonIncrSeqs(h, as, bs), function(hs) {
          c(h, hs)
        })
      }))
    }
  }
  worker <- function(rlrls, smusmus, aacc, lastx0lastrowt) {
    if(length(rlrls) <= 1L) {
      1L
    } else {
      x0 <- smusmus[1L] - aacc[1L]
      rl <- rlrls[1L]
      rls <- rlrls[-1L]
      smus <- smusmus[-1L]
      acc <- aacc[-1L]
      nacc <- length(acc)
      lastx0 <- lastx0lastrowt[1L]
      lastrowt <- lastx0lastrowt[-1L]
      aas <- vapply(c(max(lastx0, x0), lastrowt), function(i) {
        max(rl, i)
      }, integer(1L))
      rows <- boundedNonIncrSeqs(x0, aas, lambda)
      sum(vapply(rows, function(row) {
        l <- length(row) - 1L
        trow <- tail(row, l)
        irow <- head(row, l)
        m <- min(nacc, l)
        worker(rls, smus, head(acc, m) + head(trow, m), irow)
      }, integer(1L)))
    }
  }
  worker(revlam, cumsum(mu), rep(0L, n-1L), rep(0L, n))
}

#' @importFrom utils head tail
#' @noRd
.pieriRule <- function(lambda, n) {
  if(n == 0L) {
    return(list(lambda))
  }
  if(n < 0L) {
    return(list())
  }
  go <- function(k, d_ds, p_ps, l_ls) {
    elld_ds <- length(d_ds)
    if(elld_ds == 0L) {
      if(k == 0L) {
        return(list(integer(0L)))
      } else {
        return(list())
      }
    }
    d <- d_ds[1L]
    l <- l_ls[1L]
    if(elld_ds == 1L) {
      if(k <= d) {
        if((lpk <- l + k) > 0L) {
          return(list(lpk))
        } else {
          return(list(integer(0L)))
        }
      } else {
        return(list())
      }
    }
    if(k > p_ps[1L]) {
      return(list())
    }
    ps <- p_ps[-1L]
    q <- ps[1L]
    ds <- d_ds[-1L]
    ls <- l_ls[-1L]
    as <- .rg(max(0L, k-q), min(d, k))
    do.call(
      c,
      lapply(as, function(a) {
        h <- l + a
        tls <- go(k-a, ds, ps, ls)
        lapply(tls, function(tl) {
          c(h, tl)
        })
      })
    )
  }
  ellLambda <- length(lambda)
  if(ellLambda == 0L) {
    diffs <- n
  } else if(ellLambda == 1L) {
    diffs <- c(n, lambda)
  } else {
    diffs <- c(n, head(lambda, -1L) - tail(lambda, -1L), lambda[ellLambda])
  }
  dsums <- rev(cumsum(rev(diffs)))
  go(n, diffs, dsums, c(lambda, 0L))
}

.iteratedPieriRulepp <- function(plambda, coeff0, ns) {
  f <- function(t0, c_ps) {
    c <- c_ps[[1L]]
    ps <- c_ps[[2L]]
    g <- function(t, p) {
      pAsString <- partitionAsString(p)
      if(pAsString %in% names(t)) {
        t[[pAsString]][["value"]] <- t[[pAsString]][["value"]] + c
      } else {
        t[[pAsString]] <- list("lambda" = p, "value" = c)
      }
      t
    }
    Reduce(g, ps, init = t0)
  }
  worker <- function(old, n_ns) {
    if(length(n_ns) == 0L) {
      return(old)
    }
    n <- n_ns[1L]
    stuff <- lapply(old, function(lam_coeff) {
      list(lam_coeff[["value"]], .pieriRule(lam_coeff[["lambda"]], n))
    })
    new <- Reduce(f, stuff, init = list())
    worker(new, n_ns[-1L])
  }
  worker(list(list("lambda" = plambda, "value" = coeff0)), ns)
}

.iteratedPieriRulep <- function(plambda, ns) {
  .iteratedPieriRulepp(plambda, 1L, ns)
}

.iteratedPieriRule <- function(ns) {
  .iteratedPieriRulep(integer(0L), ns)
}

.KostkaNumbersWithGivenMu <- function(mu) {
  .iteratedPieriRule(rev(mu))
}

#' @title Kostka numbers with given \eqn{\mu}
#' @description Lists all positive Kostka numbers \eqn{K(\lambda,\mu)} with 
#'   a given partition \eqn{\mu}.
#' 
#' @param mu integer partition
#' @param output the format of the output, either \code{"vector"} or 
#'   \code{"list"}
#'
#' @return If \code{output="vector"}, this function returns a named vector. 
#'   This vector is made of the positive Kostka numbers 
#'   \eqn{K(\lambda,\mu)} and its names encode the partitions \eqn{\lambda}.
#'   If \code{ouput="list"}, this function returns a list of lists. 
#'   Each of these lists has two 
#'   elements. The first one is named \code{lambda} and is an integer 
#'   partition, and the second one is named \code{value} and is a positive 
#'   integer, the Kostka number \eqn{K(\lambda,\mu)}. It is faster to 
#'   compute the Kostka numbers with this function than computing the 
#'   individual Kostka numbers with the function \code{\link{KostkaNumber}}.
#' @export
#' @seealso \code{\link{KostkaNumber}}, 
#'   \code{\link{KostkaNumbersWithGivenLambda}}.
#'
#' @examples
#' KostkaNumbersWithGivenMu(c(2, 1, 1))
KostkaNumbersWithGivenMu <- function(mu, output = "vector") {
  stopifnot(isPartition(mu))
  output <- match.arg(output, c("vector", "list"))
  kNumbers <- .KostkaNumbersWithGivenMu(removeTrailingZeros(as.integer(mu)))
  if(output == "vector") {
    vapply(kNumbers, `[[`, integer(1L), "value", USE.NAMES = TRUE)
  } else {
    kNumbers
  }
}

#' @title Kostka numbers with given \eqn{\lambda}
#' @description Lists all positive Kostka numbers \eqn{K(\lambda,\mu)} with 
#'   a given partition \eqn{\lambda}.
#' 
#' @param lambda integer partition
#' @param output the format of the output, either \code{"vector"} or 
#'   \code{"list"}
#'
#' @return If \code{output="vector"}, this function returns a named vector. 
#'   This vector is made of the non-zero (i.e. positive) Kostka numbers 
#'   \eqn{K(\lambda,\mu)}, which are integers, and its names encode the 
#'   partitions \eqn{\mu}.
#'   If \code{ouput="list"}, this function returns a list of lists. 
#'   Each of these lists has two 
#'   elements. The first one is named \code{mu} and is an integer 
#'   partition, and the second one is named \code{value} and is a positive 
#'   integer, the Kostka number \eqn{K(\lambda,\mu)}. It is faster to 
#'   compute the Kostka numbers with this function than computing the 
#'   individual Kostka numbers with the function \code{\link{KostkaNumber}}.
#' @export
#' @importFrom partitions conjugate
#' @importFrom utils tail
#' @seealso \code{\link{KostkaNumber}}, 
#'   \code{\link{KostkaNumbersWithGivenMu}}.
#'
#' @examples
#' KostkaNumbersWithGivenLambda(c(2, 1, 1))
KostkaNumbersWithGivenLambda <- function(lambda, output = "vector") {
  stopifnot(isPartition(lambda))
  output <- match.arg(output, c("vector", "list"))
  lambda <- removeTrailingZeros(as.integer(lambda))
  mus <- rev(.dominatedPartitions(lambda))
  nmus <- length(mus)
  kNumbers <- rep(1L, nmus)
  musAsStrings <-
    vapply(mus, partitionAsString, character(1L), USE.NAMES = FALSE)
  names(kNumbers) <- musAsStrings
  ellLambda <- length(lambda)
  if(nmus >= 2L && ellLambda >= 2L) {
    names(mus) <- musAsStrings
    lambdap <- conjugate(lambda)
    nlambda <- sum(seq_len(ellLambda - 1L) * tail(lambda, -1L))
    nlambdap <- sum(seq_len(lambda[1L] - 1L) * tail(lambdap, -1L))
    elambda <- nlambdap - nlambda
    for(muAsString in tail(musAsStrings, -1L)) {
      mu <- mus[[muAsString]]
      mup <- conjugate(mu)
      ellMu <- mup[1L]
      i_ <- seq_len(ellMu - 1L)
      nmu <- sum(i_ * tail(mu, -1L))
      nmup <- sum(seq_len(mu[1L] - 1L) * tail(mup, -1L))
      emu <- nmup - nmu
      ee <- elambda - emu 
      x <- 0L
      for(i in i_) {
        mu_i <- mu[i]
        for(j in (i+1L):ellMu) {
          mu_j <- mu[j]
          dmuij <- mu_i - mu_j
          kappa <- mu
          for(t in seq_len(mu_j - 1L)) {
            kappa[i] <- mu_i + t
            kappa[j] <- mu_j - t
            kappaOrd <- sort(kappa, decreasing = TRUE)
            kappaOrdAsString <- partitionAsString(kappaOrd)
            if(kappaOrdAsString %in% musAsStrings){
              x <- x + kNumbers[kappaOrdAsString] * (dmuij + 2L*t) 
            }
          }
          mu_i_plus_mu_j <- mu_i + mu_j
          kappa[i] <- mu_i_plus_mu_j
          kappaOrd <- sort(kappa[-j], decreasing = TRUE)
          kappaOrdAsString <- partitionAsString(kappaOrd)
          if(kappaOrdAsString %in% musAsStrings){
            x <- x + kNumbers[kappaOrdAsString] * mu_i_plus_mu_j
          }
        }
      }
      kNumbers[muAsString] <- x %/% ee
    }
  }
  if(output == "list") {
    kNumbers <- mapply(
      function(kNumber, mu) {
        list("mu" = mu, "value" = kNumber)
      },
      kNumbers, mus,
      USE.NAMES = TRUE, SIMPLIFY = FALSE
    )
  }
  kNumbers
}

#' @title Kostka numbers for all partitions of a given weight
#' @description Computes the Kostka numbers for all integer partitions 
#'   of a given weight
#' @param n positive integer, the weight of the partitions
#'
#' @return An integer matrix, whose row names and column names encode the 
#'   partitions \eqn{\lambda} and \eqn{\mu} and whose entries are the 
#'   Kostka numbers \eqn{K(\lambda,\mu)}.
#' @export
#'
#' @examples
#' KostkaNumbers(4)
KostkaNumbers <- function(n) {
  stopifnot(isPositiveInteger(n))
  if(n == 0L) {
    Knumbers <- matrix(1L)
    colnames(Knumbers) <- rownames(Knumbers) <- partitionAsString(integer(0L))
    return(Knumbsers)
  }
  lambdas <- listOfPartitions(n)
  lambdasAsStrings <-
    vapply(lambdas, partitionAsString, character(1L), USE.NAMES = FALSE)
  zeros <- rep(0L, length(lambdas))
  names(zeros) <- lambdasAsStrings
  Knumbers <- do.call(
    rbind,
    lapply(lambdas, function(lambda) {
      kNumbersLambda <-
        KostkaNumbersWithGivenLambda(lambda, output = "vector")
      zeros[names(kNumbersLambda)] <- kNumbersLambda
      zeros
    })
  )
  colnames(Knumbers) <- rownames(Knumbers) <- lambdasAsStrings
  Knumbers
}