isTableau <- function(tableau) {
  if(!is.list(tableau) || !isPartition(lengths(tableau))) {
    return(FALSE)
  }
  modes <- vapply(tableau, mode, character(1L))
  length(unique(modes)) == 1L
}

modeTableau <- function(tableau) {
  unique(vapply(tableau, mode, character(1L)))
}

#' @title Shape of a tableau
#' @description The shape of a tableau.
#'
#' @param tableau a tableau (list of vectors having the same mode)
#'
#' @return The shape of the tableau (an integer partition).
#' @export
#'
#' @examples
#' tableau <- list(c(2, 1, 3), c(5, 2))
#' tableauShape(tableau)
tableauShape <- function(tableau) {
  stopifnot(isTableau(tableau))
  lengths(tableau)
}

.dualTableau <- function(tableau) {
  lapply(seq_along(tableau[[1L]]), function(i) {
    row <- sapply(tableau, `[`, i)
    row[!is.na(row)]
  })
}

#' @title Dual tableau
#' @description The dual tableau of a tableau (mirror image to the main 
#'   diagonal).
#' 
#' @param tableau a tableau
#'
#' @return A tableau.
#' @export
#'
#' @examples
#' tbl <- list(c("a", "s", "e", "f"), c("f", "o"), c("u"))
#' dualTableau(tbl)
dualTableau <- function(tableau) {
  stopifnot(isTableau(tableau))
  .dualTableau(tableau)
}

#' @title Pretty tableau
#' @description Pretty form of a tableau.
#'
#' @param tableau a tableau
#'
#' @return A '\code{noquote}' character matrix.
#' @export
#'
#' @examples
#' tbl <- list(c(0, 2, 1, 1), c(4, 1), c(1, 2))
#' prettyTableau(tbl)
prettyTableau <- function(tableau) {
  stopifnot(isTableau(tableau))
  ls <- lengths(tableau)
  n <- ls[1L]
  M <- t(vapply(seq_along(ls), function(i) {
    c(as.character(tableau[[i]]), rep("", n - ls[i]))
  }, character(n)))
  rownames(M) <- paste0(seq_len(nrow(M)), " ->")
  colnames(M) <- rep("", ncol(M))
  noquote(M)
}
