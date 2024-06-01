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
      do.call(c, lapply(rows, function(row) {
        worker(
          rls, 
          smus, 
          acc + head(tail(row, -1L), -1L), 
          head(row, -1L), 
          c(list(row), table)
        )
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
    } else {
      list(integer(0L))
    }
  }
  worker(revLambda, cumsum(mu), rep(0L, n-1L), rep(0L, n), list())
}

skewPartitionRows <- function(lambda, mu) {
  ellLambda <- length(lambda)
  ellMu <- length(mu)
  mu <- c(mu, rep(0L, ellLambda - ellMu))
  unlist(lapply(seq_len(ellLambda), function(i) {
    rep(i, lambda[i] - mu[i])
  }))  
}

# mkSkewPartition <- function(lambda, mu) {
#   ellLambda <- length(lambda)
#   ellMu <- length(mu)
#   mu <- c(mu, rep(0L, ellLambda - ellMu))
#   cbind(mu, lambda - mu)
# }

.adjust <- function(f, i, tableau) {
  row <- tableau[[i]]
  tableau[[i]] <- f(row)
  tableau
}

.growTableau <- function(j, tableau, lambda, mu) {
  f <- function(tbl, i) {
    .adjust(
      function(row) {
        c(row, j)
      },
      i,
      tbl
    )
  }
  Reduce(f, as.list(skewPartitionRows(lambda, mu)), init = tableau, simplify = FALSE)
}

GTpatternToTableau <- function(pattern) {
  l <- length(pattern)
  diagonals <- lapply(seq_len(l), function(j) {
    removezeros(vapply(seq_len(j), function(i) {
      pattern[[l-i+1L]][i]
    }, integer(1L)))
  })
  lambda <- diagonals[[l]]
  ellLambda <- length(lambda)
  startingTableau <- replicate(m, integer(0L), simplify = FALSE)
  go <- function(i, tableau) {
    if(i == 0L) {
      go(
        1L,
        .adjust(
          function(row) {
            rep(1L, diagonals[[1L]])
          },
          1L,
          tableau
        )
      )
    } else if(i < l) {
      go(
        i + 1L,
        .growTableau(i + 1L, tableau, diagonals[[i+1L]], diagonals[[i]])
      )
    } else {
      tableau
    }
  }
  go(0L, startingTableau)
}

# gtPatternToTableau pattern = 
#   if l >= 0 
#     then DF.toList $ go 0 startingTableau
#     else [S.singleton 1]
#   where
#     (corner, diagonals) = gtPatternDiagonals pattern
#     diagonals' = toPartitionUnsafe [corner] : diagonals
#     l = length diagonals - 1
#     lambda = diagonals !! l
#     m = partitionWidth lambda
#     startingTableau = S.replicate m S.Empty
#     zippedDiagonals = zip diagonals diagonals'
#     skewPartition i = mkSkewPartition (zippedDiagonals !! i)
#     go i tableau
#       | i == 0 = go 1 (S.adjust' (flip (><) (S.replicate corner 1)) 0 tableau)
#       | i == l+2 = tableau
#       | otherwise = 
#           go (i+1) (growTableau (i+1) tableau (skewPartition (i-1)))
#     growTableau :: Int -> Seq (Seq Int) -> SkewPartition -> Seq (Seq Int)
#     growTableau j tableau skewPart =
#       DF.foldr (\(i, _) -> S.adjust' (flip (|>) j) (i-1)) tableau 
#                 (skewPartitionElements skewPart)
