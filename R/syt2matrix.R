.tableau2matrix <- function(tableau){
  Matrix::sparseMatrix(i = rep(seq_along(tableau), times = lengths(tableau)),
                       j = unlist(sapply(lengths(tableau), seq_len)),
                       x = unlist(tableau))
}

#' Tableau as sparse matrix
#' @description Representation of a tableau as a sparse matrix; only for a 
#'   tableau with numeric or logical entries.
#'
#' @param tableau a tableau with numeric or logical entries
#'
#' @return A sparse matrix.
#' @importFrom Matrix sparseMatrix
#' @export
#'
#' @examples
#' syt <- list(c(1, 2, 6), c(3, 5), 4)
#' tableau2matrix(syt)
tableau2matrix <- function(tableau){
  stopifnot(isTableau(tableau))
  if(!is.element(modeTableau(tableau), c("numeric", "logical"))) {
    stop(
      "This function only applies to tableaux with numeric or logical entries."
    )
  }
  .tableau2matrix(tableau)
}


#' Standard Young tableau as sparse matrix
#' @description Representation of a standard Young tableau as a sparse matrix.
#'
#' @param syt a standard Young tableau
#'
#' @return A sparse matrix.
#' @seealso \code{\link{matrix2syt}}.
#' @importFrom Matrix sparseMatrix
#' @export
#' 
#' @note This function is the same as \code{\link{tableau2matrix}} except that 
#'   in addition it checks that the given tableau is a standard Young tableau.
#'
#' @examples
#' syt <- list(c(1, 2, 6), c(3, 5), 4)
#' syt2matrix(syt)
syt2matrix <- function(syt){
  checkSYT(syt)
  .tableau2matrix(syt)
}


.matrix2syt <- function(M) {
  sapply(seq_len(nrow(M)),
         function(i) as.integer(removezeros(M[i, ])), simplify = FALSE)
}

#' Standard Young tableau from a matrix
#' @description Converts a matrix to a standard Young tableau.
#'
#' @param M a matrix
#'
#' @return A standard Young tableau.
#' @seealso \code{\link{syt2matrix}}.
#' @export
#'
#' @examples
#' M <- rbind(c(1,2,6), c(3,5,0), c(4,0,0))
#' matrix2syt(M)
matrix2syt <- function(M){
  out <- .matrix2syt(M)
  if(isSYT(out)){
    return(out)
  }else{
    stop("Invalid matrix.")
  }
}
