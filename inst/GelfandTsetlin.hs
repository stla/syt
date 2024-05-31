-- | Generates all Kostka-Gelfand-Tsetlin tableau, that is, triangular arrays like
--
-- > [ 3 ]
-- > [ 3 , 2 ]
-- > [ 3 , 1 , 0 ]
-- > [ 2 , 0 , 0 , 0 ]
--
-- with both rows and column non-increasing such that
-- the top diagonal read lambda (in this case @lambda=[3,2]@) and the diagonal sums
-- are partial sums of mu (in this case @mu=[2,1,1,1]@)
--
-- The number of such GT tableaux is the Kostka
-- number K(lambda,mu).
--
kostkaGelfandTsetlinPatterns' :: Partition -> [Int] -> [GT]
kostkaGelfandTsetlinPatterns' plam@(Partition lambda0) mu0
  | minimum mu0 < 0                       = []
  | wlam == 0                             = if wmu == 0 then [ [] ] else []
  | wmu  == wlam && plam `dominates` pmu  = list
  | otherwise                             = []
  where

    pmu = mkPartition mu0

    nlam = length lambda0
    nmu  = length mu0

    n = max nlam nmu

    lambda = lambda0 ++ replicate (n - nlam) 0
    mu     = mu0     ++ replicate (n - nmu ) 0

    revlam = reverse lambda

    wmu  = sum' mu
    wlam = sum' lambda

    list = worker 
             revlam 
             (scanl1 (+) mu) 
             (replicate (n-1) 0) 
             (replicate (n  ) 0) 
             []

    worker
      :: [Int]       -- lambda_i in reverse order
      -> [Int]       -- partial sums of mu
      -> [Int]       -- sums of the tails of previous rows
      -> [Int]       -- last row
      -> [[Int]]     -- the lower part of GT tableau we accumulated so far (this is not needed if we only want to count)
      -> [GT]   

    worker (rl:rls) (smu:smus) (a:acc) (lastx0:lastrowt) table = stuff 
      where
        x0 = smu - a
        stuff = concat 
          [ worker rls smus (zipWith (+) acc (tail row)) (init row) (row:table)
          | row <- boundedNonIncrSeqs' x0 (map (max rl) (max lastx0 x0 : lastrowt)) lambda
          ]
    worker [rl] _ _ _ table = [ [rl]:table ] 
    worker []   _ _ _ _     = [ []         ]

    boundedNonIncrSeqs' :: Int -> [Int] -> [Int] -> [[Int]]
    boundedNonIncrSeqs' = go where
      go h0 (a:as) (b:bs) = [ h:hs | h <- [(max 0 a)..(min h0 b)] , hs <- go h as bs ]
      go _  []     _      = [[]]
      go _  _      []     = [[]]
      