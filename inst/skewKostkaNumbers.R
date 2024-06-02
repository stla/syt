skewKostkaNumbers <- function(lambda, mu, output = "vector") {
  output <- match.arg(output, c("vector", "list"))
  partitions <- 
    apply(parts(sum(lambda) - sum(mu)), 2L, syt:::removezeros, simplify = FALSE)
  names(partitions) <- 
    vapply(partitions, qspray:::partitionAsString, character(1L))
  lr <- jack::LRskew(lambda, mu, output = "list")
  nus <- lr[["nu"]]
  coeffs <- lr[["coeff"]]
  kNumbers <- vapply(partitions, function(partition) {
    sum(coeffs * vapply(nus, function(nu) {
      KostkaNumber(nu, partition)
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

