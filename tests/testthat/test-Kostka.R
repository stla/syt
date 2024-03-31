test_that("Kostka numbers", {
  lambda <- c(4, 2, 1); n <- 5
  nssytx <- count_ssytx(lambda, n)
  compos <- partitions::compositions(sum(lambda), n)
  expect_true(
    nssytx == sum(apply(compos, 2L, function(w) KostkaNumber(lambda, w)))
  )
})