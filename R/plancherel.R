connectedParts <- function(lambda0){
  lambda <- lambda0
  lambda[1L] <- lambda[1L]+1L
  out <- list(lambda)
  for(i in seq_along(lambda0)[-1L]){
    if(lambda0[i]<lambda0[i-1L]){
      lambda <- lambda0
      lambda[i] <- lambda[i]+1L
      out <- c(out, list(lambda))
    }
  }
  c(out, list(c(lambda0,1L)))
}

#' Plancherel growth process
#' @description Samples a path of the Young graph according to the Plancherel
#' growth process.
#'
#' @param n the size of the path to be sampled
#'
#' @return The path as a list, starting from the root vertex \code{1}.
#' @export
#'
#' @seealso \code{\link{gprocess2syt}} and \code{\link{syt2gprocess}} to convert
#' a Young path to a standard Young tableau and conversely.
#'
#' @examples
#' rgprocess(7)
rgprocess <- function(n){
  path <- list(1L)
  for(i in seq_len(n-1L)){
    set <- connectedParts(path[[i]])
    f_nu <- count_sytx(path[[i]])
    probs <- sapply(set, function(p) count_sytx(p))/(i+1)/f_nu
    path <- c(path, sample(set, 1L, prob=probs))
  }
  path
}

