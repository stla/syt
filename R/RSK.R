bump <- function(P, Q, e, i){
  if(length(P)==0L) return(list(P=list(e), Q=list(i)))
  p <- P[[1L]]
  if(e > p[length(p)]){
    P[[1L]] <- c(p, e); Q[[1L]] <- c(Q[[1L]], i)
    return(list(P=P, Q=Q))
  }else{
    j <- which.min(p < e)
    w <- p[j]; P[[1L]][j] <- e
    b <- bump(P[-1L], Q[-1L], w, i)
    return(list(P=c(P[1L], b$P), Q=c(Q[1L], b$Q)))
  }
}

#' Robinson-Schensted correspondence
#' @description Pair of standard Young tableaux given from a permutation
#' by the Robinson-Schensted correspondence.
#'
#' @param sigma a permutation given as a vector of integers
#'
#' @return A list of two standard Young tableaux.
#' @export
#'
#' @examples
#' RS(c(1, 3, 6, 4, 7, 5, 2))
RS <- function(sigma){
  if(!setequal(sigma,seq_along(sigma))){
    stop("`sigma` is not a permutation", call. = FALSE)
  }
  sigma <- as.integer(sigma)
  out <- bump(list(), list(), sigma[1L], 1L)
  for(i in seq_along(sigma)[-1L]){
    out <- bump(out$P, out$Q, sigma[i], i)
  }
  return(out)
}
