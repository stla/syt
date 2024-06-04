xxx <- function(lambda, mu, w) {
  # if(length(w) == 0L) {
  #   return(list(list(mu)))
  # }
  # if(sum(lambda) == sum(mu)) {
  #   
  # }
  d <- sum(lambda) - w[length(w)]
  if(d == sum(mu)) {
    mu <- c(mu, rep(0L, length(lambda)-length(mu)))
    if(all(lambda >= mu) &&
       all(head(mu, -1L) >= tail(lambda, -1L))
    ) {
      return(list(list(lambda,mu)))
    } else {
      return(list())
    }
  }
  ps <- syt:::.partitionsFittingRectangleWithZeros(
    lambda[1], length(lambda), d
  )
  ps <- Filter(
    function(kappa) {
      all(lambda >= kappa) &&
        all(head(kappa, -1L) >= tail(lambda, -1L))
    }, 
    ps
  )
  do.call(
    c,
    lapply(ps, function(kappa) {
      lapply(xxx(kappa, mu, head(w, -1)), function(x) {
        c(list(lambda),x)
      })
    })  
  )
}

lambda <- c(4L, 3L, 3L, 2L, 1L, 1L)

boundedNonIncrSeqs <- function(h0, a_as, b_bs) {
  laas <- length(a_as)
  lbbs <- length(b_bs)
  if(laas >= 1L && lbbs >= 1L) {
    a <- a_as[1L]
    as <- a_as[-1L]
    b <- b_bs[1L]
    bs <- b_bs[-1L]
    hrange <- syt:::.rg(max(0L, a), min(h0, b))
    do.call(c, lapply(hrange, function(h) {
      lapply(boundedNonIncrSeqs(h, as, bs), function(hs) {
        c(h, hs)
      })
    }))
  } else {
    list(integer(0L))
  }
}

sandwichedPartitions <- function(d, h0, a_as, b_bs) {
  laas <- length(a_as)
  lbbs <- length(b_bs)
  if(d < 0L) {
    list()
  } else if(laas >= 1L && lbbs >= 1L) {
    if(d < sum(a_as) || d > sum(b_bs)) {
      return(list())
    } 
    if(d == 0L) {
      return(list(rep(0L, laas)))
    }
    a <- a_as[1L]
    as <- a_as[-1L]
    b <- b_bs[1L]
    bs <- b_bs[-1L]
    hrange <- syt:::.rg(max(0L, a), min(h0, b))
    union(
      do.call(c, lapply(hrange, function(h) {
        lapply(sandwichedPartitions(d-h, h, as, bs), function(hs) {
          c(h, hs)
        })
      })),
      list()
    )
  } else {
    if(d == 0L) {
      list(integer(0L))  
    } else {
      list()
    }
  }
}
