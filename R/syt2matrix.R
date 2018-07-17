.syt2matrix <- function(syt){
  Matrix::sparseMatrix(i = rep(seq_along(syt), times=lengths(syt)),
                       j = unlist(sapply(lengths(syt), seq_len)),
                       x = as.integer(unlist(syt)))
}

#' Standard Young tableau as sparse matrix
#' @description Representation of a standard Young tableau as a sparse matrix.
#'
#' @param syt a standard Young tableau
#'
#' @return A sparse matrix.
#' @seealso \code{\link{matrix2syt}}
#' @importFrom Matrix sparseMatrix
#' @export
#'
#' @examples
#' syt <- list(c(1,2,6), c(3,5), 4)
#' syt2matrix(syt)
syt2matrix <- function(syt){
  checkSYT(syt)
  .syt2matrix(syt)
}


.matrix2syt <- function(M){
  sapply(seq_len(nrow(M)),
         function(i) as.integer(removezeros(M[i,])), simplify = FALSE)
}

#' Standard Young tableau from a matrix
#' @description Converts a matrix to a standard Young tableau.
#'
#' @param M a matrix
#'
#' @return A standard Young tableau.
#' @seealso \code{\link{syt2matrix}}
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
