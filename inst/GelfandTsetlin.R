# lambda and mu are clean
isDominatedBy <- function(mu, lambda) {
  n <- sum(lambda)
  lambda <- c(lambda, rep(0L, n - length(lambda)))
  dominated <- TRUE
  i <- 1L
  ellMu <- length(mu)
  while(dominated && i <= ellMu) {
    dominated <- sum(head(mu, i)) <= sum(head(lambda, i))
    i <- i + 1L
  }
  dominated
}

# lambda and mu are clean
.GelfandTsetlinPatterns <- function(lambda, mu) {
  ellLambda <- length(lambda)
  ellMu <- length(mu)
  if(ellLambda == 0L) {
    if(ellMu == 0L) {
      return(list(list()))
    } else {
      return(list())
    }
  }
  wLambda <- sum(lambda)
  wMu <- sum(mu)
  if(wMu != wLambda || !isDominatedBy(mu, lambda)) {
    return(list())
  }
  n <- max(ellLambda, ellMu)
  if(ellLambda != n) {
    lambda <- c(lambda, rep(0L, n - ellLambda))
  } else if(ellMu != n) {
    mu <- c(mu, rep(0L, n - ellMu))
  } 
  revLambda <- rev(lambda)
  worker <- function(
      rl_rls,
      smu_smus,
      a_acc,
      r1_rowt,
      table
  ) {
    l <- length(rl_rls) 
    if(l >= 2L) {
      rl <- rl_rls[1L]
      rls <- rl_rls[-1L]
      smu <- smu_smus[1L]
      smus <- smu_smus[-1L]
      a <- a_acc[1L]
      acc <- a_acc[-1L]
      r1 <- r1_rowt[1L]
      rowt <- r1_rowt[-1L]
      x0 <- smu - a
      rows <- boundedNonIncrSeqs(
        x0, 
        vapply(c(max(r1, x0), rowt), function(i) {
          max(rl, i)
        }, integer(1L)),
        lambda
      )
      # print(rows)
      do.call(c, lapply(rows, function(row) {
        worker(rls, smus, acc + head(tail(row, -1L), -1L), head(row, -1L), c(list(row), table))
      }))
    } else if(l == 1L) {
      list(c(list(rl_rls), table))
    } else {
      list(list())
    }
  }
  boundedNonIncrSeqs <- function(h0, a_as, b_bs) {
    laas <- length(a_as)
    lbbs <- length(b_bs)
    if(laas >= 1L && lbbs >= 1L) {
      a <- a_as[1L]
      as <- a_as[-1L]
      b <- b_bs[1L]
      bs <- b_bs[-1L]
      hrange <- .rg(max(0L, a), min(h0, b))
      do.call(c, lapply(hrange, function(h) {
        lapply(boundedNonIncrSeqs(h, as, bs), function(hs) {
          c(h, hs)
        })
      }))
      
      # mapply(
      #   function(h, hs) {
      #     c(h, hs)
      #   },
      #   hrange,
      #   lapply(hrange, function(h) {
      #     boundedNonIncrSeqs(h, as, bs)
      #   }),
      #   SIMPLIFY = FALSE, USE.NAMES = FALSE
      # )
    } else {
      list(integer(0L))
    }
  }
  worker(revLambda, cumsum(mu), rep(0L, n-1L), rep(0L, n), list())
}
