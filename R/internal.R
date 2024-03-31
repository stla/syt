#' @importFrom utils head
removezeros <- function(x){ # e.g c(3,1,0,0) -> c(3,1)
  i <- match(0L, x)
  if(!is.na(i)) head(x, i-1L) else x
}

isPositiveInteger <- function(x) {
  x >= 1 && floor(x) == x
}

arePositiveIntegers <- function(x){
  all(x > 0 & floor(x) == x)
}

areNonnegativeIntegers <- function(x){
  all(x >= 0 & floor(x) == x)
}

isPartition <- function(x){
  is.numeric(x) && areNonnegativeIntegers(x) && all(diff(x) <= 0)
}

checkPartition <- function(x){
  if(isPartition(x)){
    return(as.integer(removezeros(x)))
  }else{
    stop("`lambda` is not a partition.", call. = FALSE)
  }
}

isStrictlyIncreasing <- function(x) {
  all(diff(x) > 0)
}

isWeaklyIncreasing <- function(x) {
  all(diff(x) >= 0)
}

checkSYTrows <- function(syt){
  all(vapply(syt, isStrictlyIncreasing, logical(1L)))
}

checkSSYTrows <- function(ssyt){
  all(vapply(ssyt, isWeaklyIncreasing, logical(1L)))
}

checkSYT <- function(syt){
  if(!isSYT(syt)){
    stop("Not a standard Young tableau.", call. = FALSE)
  }else{
    return(invisible())
  }
}

.rg <- function(start, end) {
  #seq_len(end - start + 1L) + (start - 1L)
  if(start <= end) start:end else integer(0L)
}
