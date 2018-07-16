.dsyt <- function(syt){
  .matrix2syt(t(as.matrix(.syt2matrix(syt))))
}

#' Dual tableau
#' @description The dual standard Young tableau of a standard Young tableau.
#'
#' @param syt standard Young tableau
#'
#' @return A standard Young tableau.
#' @export
#'
#' @examples
#' syt <- list(c(1,2,6), c(3,5), 4)
#' dsyt(syt)
dsyt <- function(syt){
  lapply(matrix2syt(t(as.matrix(syt2matrix(syt)))), as.integer)
}
