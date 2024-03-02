#' @title Checks whether a tableau is standard
#' @description Checks whether a tableau is a standard Young tableau.
#' 
#' @param tableau a tableau
#'
#' @return A Boolean value.
#' @export
#'
#' @examples
#' tbl <- list(c(1,2,6), c(3,5), 4)
#' isSYT(tbl)
isSYT <- function(tableau){
  contents <- unlist(tableau)
  N <- length(contents)
  isTableau(tableau) &&
    modeTableau(tableau) == "numeric" &&
    setequal(contents, 1L:N) &&
    checkSYTrows(tableau) &&
    checkSYTrows(.dualTableau(tableau))
}

#' @title Checks whether a tableau is semistandard
#' @description Checks whether a tableau is a semistandard Young tableau.
#' 
#' @param tableau a tableau
#'
#' @return A Boolean value.
#' @export
#'
#' @examples
#' tbl <- list(c(1,2,6), c(5,5), 7)
#' isSSYT(tbl)
isSSYT <- function(tableau){
  stopifnot(isTableau(tableau))
  contents <- unlist(tableau)
  arePositiveIntegers(contents - 1L) &&
    checkSSYTrows(tableau) && 
    checkSYTrows(.dualTableau(tableau))
}
