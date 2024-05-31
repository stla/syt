#' Kostka numbers
#' @description Computes a Kostka number.
#'
#' @param lambda an integer partition
#' @param mu an integer vector whose sum equals the weight (i.e. the sum) of 
#'   \code{lambda}
#'
#' @return The Kostka number corresponding to \code{lambda} and \code{mu}.
#' @export
#' @importFrom utils head tail
#'
#' @details
#' The Kostka number \eqn{K(\lambda,\mu)} is the number of semistandard 
#'   Young tableaux with shape \eqn{\lambda} and content \eqn{\mu}. It does 
#'   not depend on the order of the elements of \eqn{\mu} (so one can always 
#'   take an integer partition for \eqn{\mu}). The \emph{content} is the 
#'   vector whose \eqn{i}-th element is the number of occurrences of \eqn{i} 
#'   in the tableau.
#'   
#' @examples
#' KostkaNumber(c(3,2), c(1,1,1,2))
#' KostkaNumber(c(3,2), c(1,1,2,1))
#' KostkaNumber(c(3,2), c(1,2,1,1))
#' KostkaNumber(c(3,2), c(2,1,1,1))
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
