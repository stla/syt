#' @importFrom utils head
removezeros <- function(x){ # e.g c(3,1,0,0) -> c(3,1)
  i <- match(0L, x)
  if(!is.na(i)) head(x, i-1L) else x
}

isPositiveInteger <- function(x) {
  x >= 1 && floor(x) == x
}

arePositiveIntegers <- function(x){
  all(x >= 0 & floor(x) == x)
}

isPartition <- function(x){
  is.numeric(x) && arePositiveIntegers(x) && all(diff(x) <= 0)
}

checkPartition <- function(x){
  if(isPartition(x)){
    return(removezeros(x))
  }else{
    stop("`lambda` is not a partition.", call. = FALSE)
  }
}

checkSYTrows <- function(syt){
  all(sapply(syt, function(x) all(diff(x) > 0)))
}

checkSSYTrows <- function(ssyt){
  all(sapply(ssyt, function(x) all(diff(x) >= 0)))
}

checkSYT <- function(syt){
  if(!isSYT(syt)){
    stop("Not a standard Young tableau.", call. = FALSE)
  }else{
    return(invisible())
  }
}
