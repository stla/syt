#' @importFrom utils head
removezeros <- function(x){ # e.g c(3,1,0,0) -> c(3,1)
  i <- match(0L, x)
  if(!is.na(i)) head(x, i-1L) else x
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
  all(sapply(syt, function(x) all(diff(x)>0)))
}

isSYT <- function(syt){
  contents <- unlist(syt)
  N <- length(contents)
  is.list(syt) &&
    isPartition(lengths(syt)) &&
    setequal(contents, 1L:N) &&
    checkSYTrows(syt) &&
    checkSYTrows(.dualsyt(syt))
}

checkSYT <- function(syt){
  if(!isSYT(syt)){
    stop("Not a standard Young tableau.", call. = FALSE)
  }else{
    return(invisible())
  }
}
