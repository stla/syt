test_that("Kostka numbers", {
  lambda <- c(4, 2, 1); n <- 5
  nssytx <- count_ssytx(lambda, n)
  compos <- partitions::compositions(sum(lambda), n)
  expect_true(
    nssytx == sum(apply(compos, 2L, function(w) KostkaNumber(lambda, w)))
  )
})

test_that("Kostka numbers with given mu", {
  mu <- c(3, 3, 2, 1)
  kNumbers_mu <- KostkaNumbersWithGivenMu(mu, output = "vector")
  lambdas <- partitions::parts(sum(mu))
  kNumbers <- apply(lambdas, 2L, function(lambda) KostkaNumber(lambda, mu))
  names(kNumbers) <- apply(lambdas, 2L, function(lambda) {
    partitionAsString(removeTrailingZeros(lambda))
  })
  expect_true(
    all(kNumbers[names(kNumbers_mu)] == kNumbers_mu)
  )
  expect_true(
    all(kNumbers[setdiff(names(kNumbers), names(kNumbers_mu))] == 0L)
  )
})

test_that("Kostka numbers with given lambda", {
  lambda <- c(3, 3, 2, 1)
  kNumbers_lambda <- KostkaNumbersWithGivenLambda(lambda, output = "vector")
  mus <- partitions::parts(sum(lambda))
  kNumbers <- apply(mus, 2L, function(mu) KostkaNumber(lambda, mu))
  names(kNumbers) <- apply(mus, 2L, function(mu) {
    partitionAsString(removeTrailingZeros(mu))
  })
  expect_true(
    all(kNumbers[names(kNumbers_lambda)] == kNumbers_lambda)
  )
  expect_true(
    all(kNumbers[setdiff(names(kNumbers), names(kNumbers_lambda))] == 0L)
  )
})
