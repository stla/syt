checkBallot <- function(a){
  if(!arePositiveIntegers(a-1)){
    stop("Not a ballot sequence.", call.=FALSE)
  }
  for(k in seq_along(a)){
    init <- head(a,k)
    for(r in seq_len(max(init))){
      test <- length(which(init==r)) >= length(which(init==r+1L))
      if(!test){
        stop("Not a ballot sequence.", call.=FALSE)
      }
    }
  }
  return(invisible())
}

.syt2ballot <- function(syt){
  N <- sum(lengths(syt))
  a <- integer(0L)
  for(i in 1L:N){
    k <- 0L
    j <- integer(0L)
    while(length(j)==0L){
      k <- k+1L
      j <- which(syt[[k]]==i)
      if(length(j)) a <- c(a,k)
    }
  }
  a
}

.ballot2syt <- function(a){
  lapply(1L:max(a), function(i) which(a == i))
  # syt <- vector("list", max(a))
  # for(i in seq_along(a)){
  #   syt[[a[i]]] <- c(syt[[a[i]]], i)
  # }
  # syt
}

#' Tableau as ballot sequence
#' @description Converts a standard Young tableau to its corresponding ballot
#' sequence.
#'
#' @param syt standard Young tableau
#'
#' @return A ballot sequence.
#' @seealso \code{\link{ballot2syt}}
#' @export
#'
#' @examples
#' syt <- list(c(1,2,6), c(3,5), 4)
#' syt2ballot(syt)
syt2ballot <- function(syt){
  checkSYT(syt)
  .syt2ballot(syt)
}

#' Tableau as ballot sequence
#' @description Converts a ballot sequence to its corresponding standard
#' Young tableau.
#'
#' @param a ballot sequence
#'
#' @return A standard Young tableau.
#' @seealso \code{\link{syt2ballot}}
#' @export
#'
#' @examples
#' a <- c(1,1,2,3,2,1)
#' ballot2syt(a)
ballot2syt <- function(a){
  checkBallot(a)
  .ballot2syt(a)
}
