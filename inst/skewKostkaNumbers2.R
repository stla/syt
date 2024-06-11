lastSubpartition <- function(w, lambda) {
  if(length(lambda) == 0L) {
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

skewKostkaNumbers2 <- function(lambda, mu, output = "vector") {
  output <- match.arg(output, c("vector", "list"))
  kappa <- lastSubpartition(sum(lambda)-sum(mu), lambda)
  nus <- apply(
    jack:::dominatedPartitions(kappa), 2L, syt:::removeTrailingZeros, simplify = FALSE
  )
  names(nus) <- 
    vapply(nus, qspray:::partitionAsString, character(1L))
  lr <- jack::LRskew(lambda, mu, output = "list")
  pis <- lr[["nu"]]
  coeffs <- lr[["coeff"]]
  dpis <- lapply(pis, dominatedPartitions)
  kNumbers <- vapply(nus, function(nu) {
    # i_ <- vapply(pis, function(pi) {
    #   jack:::isDominated(nu, pi)
    # }, logical(1L))
    i_ <- vapply(dpis, function(dpi) {
      is.element(list(nu), dpi)
    }, logical(1L))
    sum(coeffs[i_] * vapply(pis[i_], function(pi) {
      KostkaNumber(pi, nu)
    }, integer(1L), USE.NAMES = FALSE))
  }, integer(1L), USE.NAMES = TRUE)
  if(output == "vector") {
    kNumbers
  } else {
    mapply(
      function(nu, value) {
        list("nu" = nu, "value" = value)
      },
      partitions, kNumbers,
      USE.NAMES = TRUE, SIMPLIFY = FALSE
    )
  }
}

lambda <- c(4, 3, 3, 2, 1, 1)
mu <- c(2, 2, 1)
w <- c(3, 3, 2, 1)
skewKostkaNumbers(lambda, mu)

dominatingPartitions <- function(nu) {
  ps <- apply(partitions::parts(sum(nu)), 2L, syt:::removeTrailingZeros, simplify = FALSE)
  Filter(function(p) jack:::isDominated(nu,p), ps)
}
dd <- lapply(nus, dominatingPartitions)

dominatedPartitions <- function(nu) {
  ps <- apply(partitions::parts(sum(nu)), 2L, syt:::removeTrailingZeros, simplify = FALSE)
  Filter(function(p) jack:::isDominated(p,nu), ps)
}
