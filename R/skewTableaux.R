diffSequence <- function(x) {
  c(diff(-x), x[length(x)])
}

# semiStandardSkewTableaux :: Int -> SkewPartition -> [SkewTableau Int]
# semiStandardSkewTableaux n (SkewPartition abs) = map SkewTableau stuff where
# 
#   stuff = worker as bs ds (repeat 1) 
#   (as,bs) = unzip abs
#   ds = _diffSequence as
#   
#   -- | @worker inner outerMinusInner innerdiffs lowerbound
#   worker :: [Int] -> [Int] -> [Int] -> [Int] -> [[(Int,[Int])]]
#   worker (a:as) (b:bs) (d:ds) lb = [ (a,this):rest 
#                                    | this <- row b 1 lb 
#                                    , let lb' = (replicate d 1 ++ map (+1) this) 
#                                    , rest <- worker as bs ds lb' ] 
#   worker []     _      _      _  = [ [] ]
# 
#   -- @row length minimum lowerbound@
#   row 0  _  _       = [[]]
#   row _  _  []      = []
#   row !k !m (!a:as) = [ x:xs | x <- [(max a m)..n] , xs <- row (k-1) x as ] 

all_ssSkewTableaux <- function(lambda, mu, n) {
  row <- function(k, m, aas) {
    if(k == 0L) {
      list(integer(0L))
    } else if(length(aas) == 0L) {
      list()
    } else {
      a <- aas[1L]
      as <- aas[-1L]
      maxam <- max(a, m)
      x_ <- if(maxam <= n) maxam:n else integer(0L)
      do.call(c, lapply(x_, function(x) {
        lapply(row(k - 1L, x, as), function(xs) {
          c(x, xs)
        })
      }))
    }
  }
  worker <- function(aas, bbs, dds, lb) {
    if(length(aas) == 0L) {
      list(list())
    } else {
      a <- aas[1L]
      as <- aas[-1L]
      b <- bbs[1L]
      bs <- bbs[-1L]
      d <- dds[1L]
      print(d)
      ds <- dds[-1L]
      do.call(c, lapply(row(b, 1L, lb), function(this) {
        lbprime <- c(rep(1L, d), this + 1L)
        lapply(worker(as, bs, ds, lbprime), function(rest) {
          c(list(list(a, this)), rest)
        })
      }))
    }
  }
  as <- c(mu, rep(0L, length(lambda) - length(mu)))
  bs <- lambda - c(mu, rep(0L, length(lambda) - length(mu)))
  ds <- diffSequence(as)
  worker(as, bs, ds, rep(1L, 100))#bs[1L]))
}