#' Kostka numbers
#' @description Computes a Kostka number.
#'
#' @param lambda an integer partition
#' @param mu an integer vector whose sum equals the weight (i.e. the sum) of 
#'   \code{lambda}
#'
#' @return The Kostka number corresponding to \code{lambda} and \code{mu}.
#' @export
#' @importFrom utils head tail
#'
#' @details
#' The Kostka number \eqn{K(\lambda,\mu)} is the number of semistandard 
#'   Young tableaux with shape \eqn{\lambda} and weight \eqn{\mu}. It does 
#'   not depend on the order of the elements of \eqn{\mu} (so one can always 
#'   take an integer partition for \eqn{\mu}). The \emph{weight} is the 
#'   vector whose \eqn{i}-th element is the number of occurrences of \eqn{i} 
#'   in the tableau.
#'   
#' @examples
#' KostkaNumber(c(3,2), c(1,1,1,2))
#' KostkaNumber(c(3,2), c(1,1,2,1))
#' KostkaNumber(c(3,2), c(1,2,1,1))
#' KostkaNumber(c(3,2), c(2,1,1,1))
KostkaNumber <- function(lambda, mu) {
  lambda <- checkPartition(lambda)
  mu     <- as.integer(mu)
  wmu  <- sum(mu)
  wlam <- sum(lambda)
  if(wlam == 0L) {
    return(as.integer(wmu == 0L))
  }
  if(wmu != wlam) { #|| !jack:::isDominated(mu, lambda)) {
    return(0L)
  }
  nlam <- length(lambda)
  nmu  <- length(mu)
  n <- max(nlam, nmu)
  lambda <- c(lambda, rep(0L, n - nlam))
  mu     <- c(mu, rep(0L, n - nmu))
  revlam <- rev(lambda)
  boundedNonIncrSeqs <- function(h0, aas, bbs) {
    if(length(aas) == 0L || length(bbs) == 0L) {
      list(integer(0L))
    } else {
      a <- aas[1L]
      b <- bbs[1L]
      as <- aas[-1L]
      bs <- bbs[-1L]
      h_ <- .rg(max(0L, a), min(h0, b))
      do.call(c, lapply(h_, function(h) {
        lapply(boundedNonIncrSeqs(h, as, bs), function(hs) {
          c(h, hs)
        })
      }))
    }
  }
  worker <- function(rlrls, smusmus, aacc, lastx0lastrowt) {
    if(length(rlrls) <= 1L) {
      1L
    } else {
      x0 <- smusmus[1L] - aacc[1L]
      rl <- rlrls[1L]
      rls <- rlrls[-1L]
      smus <- smusmus[-1L]
      acc <- aacc[-1L]
      nacc <- length(acc)
      lastx0 <- lastx0lastrowt[1L]
      lastrowt <- lastx0lastrowt[-1L]
      aas <- vapply(c(max(lastx0, x0), lastrowt), function(i) {
        max(rl, i)
      }, integer(1L))
      rows <- boundedNonIncrSeqs(x0, aas, lambda)
      sum(vapply(rows, function(row) {
        l <- length(row) - 1L
        trow <- tail(row, l)
        irow <- head(row, l)
        m <- min(nacc, l)
        worker(rls, smus, head(acc, m) + head(trow, m), irow)
      }, integer(1L)))
    }
  }
  worker(revlam, cumsum(mu), rep(0L, n-1L), rep(0L, n))
}

.pieriRule <- function(lambda, n) {
  if(n == 0L) {
    return(list(lambda))
  }
  if(n < 0L) {
    return(list())
  }
  go <- function(k, d_ds, p_ps, l_ls) {
    elld_ds <- length(d_ds)
    if(elld_ds == 0L) {
      if(k == 0L) {
        return(list(integer(0L)))
      } else {
        return(list())
      }
    }
    d <- d_ds[1L]
    l <- l_ls[1L]
    if(elld_ds == 1L) {
      if(k <= d) {
        if((lpk <- l + k) > 0L) {
          return(list(lpk))
        } else {
          return(list(integer(0L)))
        }
      } else {
        return(list())
      }
    }
    if(k > p_ps[1L]) {
      return(list())
    }
    ps <- p_ps[-1L]
    q <- ps[1L]
    ds <- d_ds[-1L]
    ls <- l_ls[-1L]
    as <- .rg(max(0L, k-q), min(d, k))
    do.call(
      c,
      lapply(as, function(a) {
        h <- l + a
        tls <- go(k-a, ds, ps, ls)
        lapply(tls, function(tl) {
          c(h, tl)
        })
      })
    )
  }
  ellLambda <- length(lambda)
  if(ellLambda == 0L) {
    diffs <- n
  } else if(ellLambda == 1L) {
    diffs <- c(n, lambda)
  } else {
    diffs <- c(n, head(lambda, -1L) - tail(lambda, -1L), lambda[ellLambda])
  }
  dsums <- rev(cumsum(rev(diffs)))
  go(n, diffs, dsums, c(lambda, 0L))
}
# -- | We assume here that @lambda@ is a partition (non-increasing sequence of /positive/ integers)! 
# _pieriRule :: [Int] -> Int -> [[Int]] 
# _pieriRule lambda n
#     | n == 0     = [lambda]
#     | n <  0     = [] 
#     | otherwise  = go n diffs dsums (lambda++[0]) 
#     where
#       diffs = n : _diffSequence lambda                 -- maximum we can add to a given row
#       dsums = reverse $ scanl1 (+) (reverse diffs)    -- partial sums of remaining total we can add
#       go !k (d:ds) (p:ps@(q:_)) (l:ls) 
#         | k > p     = []
#         | otherwise = [ h:tl | a <- [ max 0 (k-q) .. min d k ] , let h = l+a , tl <- go (k-a) ds ps ls ]
#       go !k [d]    _      [l]    = if k <= d 
#                                      then if l+k>0 then [[l+k]] else [[]]
#                                      else []
#       go !k []     _      _      = if k==0 then [[]] else []

.iteratedPieriRulepp <- function(plambda, coeff0, ns) {
  f <- function(t0, c_ps) {
    c <- c_ps[[1L]]
    ps <- c_ps[[2L]]
    g <- function(t, p) {
      pAsString <- partitionAsString(p)
      if(pAsString %in% names(t)) {
        t[[pAsString]][["coeff"]] <- t[[pAsString]][["coeff"]] + c
      } else {
        t[[pAsString]] <- list("mu" = p, "coeff" = c)
      }
      t
    }
    Reduce(g, ps, init = t0)
  }
  worker <- function(old, n_ns) {
    if(length(n_ns) == 0L) {
      return(old)
    }
    n <- n_ns[1L]
    stuff <- lapply(old, function(lam_coeff) {
      list(lam_coeff[["coeff"]], .pieriRule(lam_coeff[["mu"]], n))
    })
    new <- Reduce(f, stuff, init = list())
    worker(new, n_ns[-1L])
  }
  worker(list(list("mu" = plambda, "coeff" = coeff0)), ns)
}
# iteratedPieriRule'' :: Num coeff => (Partition,coeff) -> [Int] -> Map Partition coeff
# iteratedPieriRule'' (plambda,coeff0) ns = worker (Map.singleton plambda coeff0) ns where
#   worker old []     = old
#   worker old (n:ns) = worker new ns where
#     stuff = [ (coeff, pieriRule lam n) | (lam,coeff) <- Map.toList old ] 
#     new   = foldl' f Map.empty stuff 
#     f t0 (c,ps) = foldl' (\t p -> Map.insertWith (+) p c t) t0 ps  



