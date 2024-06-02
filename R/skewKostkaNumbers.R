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
#'
#' @examples
#' skewKostkaNumbers(c(4,2,2), c(2,2))
skewKostkaNumbers <- function(lambda, mu, output = "vector") {
  LRcoeffs <- LRskew(lambda, mu, output = "list")
  output <- match.arg(output, c("vector", "list"))
  partitions <- 
    apply(parts(sum(lambda) - sum(mu)), 2L, removezeros, simplify = FALSE)
  names(partitions) <- 
    vapply(partitions, partitionAsString, character(1L))
  nus <- LRcoeffs[["nu"]]
  coeffs <- LRcoeffs[["coeff"]]
  kNumbers <- vapply(partitions, function(partition) {
    sum(coeffs * vapply(nus, function(nu) {
      KostkaNumber(nu, partition)
    }, integer(1L), USE.NAMES = FALSE))
  }, integer(1L), USE.NAMES = TRUE)
  kNumbers <- kNumbers[kNumbers != 0L]
  if(output == "vector") {
    kNumbers
  } else {
    mapply(
      function(nu, value) {
        list("nu" = nu, "value" = value)
      },
      partitions[names(kNumbers)], kNumbers,
      USE.NAMES = TRUE, SIMPLIFY = FALSE
    )
  }
}

