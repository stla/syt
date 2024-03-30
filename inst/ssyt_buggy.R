# wrong for lambda = c(2, 2, 1) and n = 4

all_ssytx <- function(lambda, n) {
  lambda <- as.integer(checkPartition(lambda))
  n <- as.integer(n)
  stopifnot(n >= 1L)
  l <- length(lambda)
  ll <- lambda[l]
  tmp <- list(list())
  res <- list()
  while(length(tmp) > 0L) {
    T <- tmp[[1L]]
    tmp <- tmp[-1L]
    k <- length(T)
    if(k == l && length(T[[k]]) >= ll) {
      res <- c(res, list(T))
    } else {
      if(k == 0L || length(T[[k]]) %in% c(0L, lambda[k])) { # 0L ????
        if(k == 0L) {
          start <- 1L
        } else {
          start <- T[[k]][1L] + 1L
        }
        for(u in rg(start, n)) {
          U <- c(T, list(u))
          tmp <- c(tmp, list(U))
        }
      } else {
        Tk <- T[[k]]
        lk <- length(Tk)
        if(k == 1L) {
          start <- Tk[lk]
        } else {
          start <- max(Tk[lk], T[[k-1L]][lk] + 1L)
        }
        for(u in rg(start, n)) {
          U <- T
          U[[k]] <- c(Tk, u)
          tmp <- c(tmp, list(U))
        }
      }
    }
  }
  res
} 
