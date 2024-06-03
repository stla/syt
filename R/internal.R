#' @importFrom utils head
#' @noRd
removezeros <- function(x){ # e.g c(3,1,0,0) -> c(3,1)
  i <- match(0L, x)
  if(!is.na(i)) head(x, i-1L) else x
}
removeTrailingZeros <- function(x) {
  n <- length(x)
  while(x[n] == 0 && n > 0L) {
    n <- n - 1L
  }
  head(x, n)
}

isPositiveInteger <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && x >= 1 && floor(x) == x
}

isIntegerVector <- function(x) {
  is.numeric(x) && !anyNA(x) && all(x == as.integer(x))
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

# the converse of toString
fromString <- function(string) {
  as.integer(strsplit(string, ",", fixed = TRUE)[[1L]])
}

partitionAsString <- function(lambda) {
  sprintf("[%s]", toString(lambda))
}

fromPartitionAsString <- function(string) {
  string <- gsub("(\\[|\\])", "", string)
  as.integer(strsplit(string, ",", fixed = TRUE)[[1L]])
}
