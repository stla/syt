vector2syt <- function(A){
  # A[i] is the row containing i
  Y <- vector("list", max(A))
  for(i in seq_along(A)){
    Y[[A[i]]] <- c(Y[[A[i]]], i)
  }
  Y
}

syt2vector <- function(Y){
  N <- sum(lengths(Y))
  A <- integer(0L)
  for(i in 1L:N){
    k <- 0L
    j <- integer(0L)
    while(length(j)==0L){
      k <- k+1L
      j <- which(Y[[k]]==i)
      if(length(j)) A <- c(A,k)
    }
  }
  A
}

#' @importFrom utils head
removezeros <- function(x){ # e.g c(3,1,0,0) -> c(3,1)
  i <- match(0L, x)
  if(!is.na(i)) return(head(x,i-1L)) else return(x)
}

arePositiveIntegers <- function(x){
  all(x>=0 & floor(x) == x)
}

isPartition <- function(x){
  is.numeric(x) && arePositiveIntegers(x) && all(diff(x) <= 0)
}

checkPartition <- function(x){
  if(isPartition(x)){
    return(removezeros(x))
  }else{
    stop("`lambda` is not a partition", call. = FALSE)
  }
}
