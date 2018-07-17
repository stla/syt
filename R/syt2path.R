whichRow <- function(M,x){
  i <- 1L
  idx <- match(x, M[i,])
  while(is.na(idx)){
    i <- i+1L
    idx <- match(x, M[i,])
  }
  i
}

#' Tableau as growth process
#' @description Converts a standard Young tableau to its corresponding growth
#' process of partitions.
#'
#' @param syt standard Young tableau
#'
#' @return A list of integer partitions, representing a path of the Young graph
#' starting from the root vertex.
#' @seealso \code{\link{gprocess2syt}}
#' @export
#'
#' @examples
#' syt <- list(c(1,2,4), 3, 5)
#' syt2gprocess(syt)
syt2gprocess <- function(syt){
  M <- syt2matrix(syt)
  N <- sum(lengths(syt))
  path <- vector("list", N)
  path[[1L]] <- c(1L, integer(length(syt)-1L))
  for(k in seq_len(N)[-1L]){
    part <- path[[k-1L]]
    i <- whichRow(M, k)
    part[i] <- part[i]+1L
    path[[k]] <- part
  }
  lapply(path, removezeros)
}


#' Growth process to tableau
#' @description Converts a growth process of integer partitions to its
#' corresponding standard Young tableau.
#'
#' @param path a path of the Young graph from the root vertex, given as a list
#' of integer partitions
#'
#' @return A standard Young tableau.
#' @seealso \code{\link{syt2gprocess}}
#' @export
#'
#' @examples
#' path <- list(1, 2, c(2,1), c(3,1), c(3,1,1))
#' gprocess2syt(path)
gprocess2syt <- function(path){
  if(!is.list(path)){
    stop("The path must be a list of integer partitions.")
  }
  path <- lapply(path, checkPartition)
  check0 <- all(sapply(path, sum) == seq_along(path))
  if(!check0){
    stop("Invalid path. The n-th element of the path must be a partition of n.")
  }
  N <- length(path)
  if(N==1L){
    return(list(1L))
  }
  check1 <- all(diff(lengths(path)) %in% c(0,1))
  if(!check1){
    stop("Invalid path.")
  }
  l <- length(path[[N]])
  syt <- vector("list", l)
  path2 <- replicate(N, integer(l))
  for(i in 1L:N){
    path2[,i][seq_along(path[[i]])] <- path[[i]]
  }
  diffs <- diff(t(path2))
  check2 <- all(apply(diffs, 1, function(x){
    sum(x) == 1 && all(x %in% c(0,1))
  }))
  if(!check2){
    stop("Invalid path.")
  }
  syt[[1L]] <- 1L
  for(i in seq_len(N-1L)){
    idx <- match(1, diffs[i,])
    syt[[idx]] <- c(syt[[idx]], i+1L)
  }
  syt
}
