orderedMatrix <- function(M) {
  columns <- apply(M, 2L, identity, simplify = FALSE)
  o <- do.call(function(...) order(..., decreasing = TRUE), columns)
  M[o, ]
}