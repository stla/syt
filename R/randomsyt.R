#' Random standard Young tableau
#' @description Uniform sampling of a standard Young tableau of a given shape.
#'
#' @param lambda shape, an integer partition
#'
#' @return A standard Young tableau of shape \code{lambda}.
#' @export
#'
#' @examples
#' rsyt(c(7,3,1))
rsyt <- function(lambda){
  lambda <- as.integer(checkPartition(lambda))
  N <- sum(lambda)
  a <- integer(N)
  i <- k <- 0L
  while(TRUE){
    i <- i + 1L
    for(j in 1L:lambda[i]){
      a[j] <- a[j] + 1L
      k <- k + 1L
    }
    if(N <= k){
      break
    }
  }
  for(m in 1L:N){
    while(TRUE){
      i <- sample.int(a[1L], 1L)
      j <- sample.int(lambda[1L], 1L)
      if(i <= a[j] && j <= lambda[i]){
        break
      }
    }
    while(TRUE){
      ih <- a[j] + lambda[i] - i - j
      if(ih == 0L){
        break
      }
      k <- sample.int(ih, 1L)
      if(k <= lambda[i] - j){
        j <- j + k
      }
      else{
        i <- k - lambda[i] + i + j
      }
    }
    lambda[i] <- lambda[i] - 1L
    a[j] <- a[j] - 1L
    a[N-m+1L] <- i
  }
  return(.ballot2syt(a))
}

