skewKostkaNumbers <- function(lambda, mu) {
  lr <- jack::LRskew(lambda, mu, output = "list")
  nus <- lr[["nu"]]
  names(nus) <- vapply(nus, qspray:::partitionAsString, character(1L))
  coeffs <- lr[["coeff"]]
  vapply(nus, function(nu0) {
    sum(coeffs * vapply(nus, function(nu) {
      KostkaNumber(nu, nu0)
    }, integer(1L), USE.NAMES = FALSE))
  }, integer(1L), USE.NAMES = TRUE)
}

lambda <- c(4, 3, 3, 2, 1, 1)
mu <- c(2, 2, 1)
w <- c(3, 3, 2, 1)
skewKostkaNumbers(lambda, mu)

