# semiStandardYoungTableaux :: Int -> Partition -> [Tableau Int]
# semiStandardYoungTableaux n part = worker (repeat 0) shape where
#   shape = fromPartition part
#   worker _ [] = [[]] 
#   worker prevRow (s:ss) 
#     = [ (r:rs) | r <- row n s 1 prevRow, rs <- worker (map (+1) r) ss ]
# 
#   -- weekly increasing lists of length @len@, pointwise at least @xs@, 
#   -- maximum value @n@, minimum value @prev@.
#   row :: Int -> Int -> Int -> [Int] -> [[Int]]
#   row _ 0   _    _      = [[]]
#   row n len prev (x:xs) = [ (a:as) | a <- [max x prev..n] , as <- row n (len-1) a xs ]
rg <- function(a, b) if(b>=a) a:b else integer(0L)
ssytxx <- function(lambda, n) {
  row <- function(n, len, prev, xxs) {
    if(len == 0L) {
      list(integer(0L))
    } else {
      x <- xxs[[1L]]
      xs <- xxs[-1L]
      do.call(c, lapply(rg(max(x, prev), n), function(a) {
        lapply(row(n, len - 1L, a, xs), function(as) {
          c(a, as)
        })
      }))
    }
  }
  worker <- function(prevRow, sss) {
    if(length(sss) == 0L) {
      list(list())
    } else {
      s <- sss[[1L]]
      ss <- sss[-1L]
      do.call(c, lapply(row(n, s, 1L, prevRow), function(r) {
        lapply(worker(r + 1L, ss), function(rs) {
          c(list(r), rs)
        })
      }))
    }
  }
  worker(rep(0L, lambda[1L]), lambda)
}

ssytxx(c(2,2,1), 4) |> length()
