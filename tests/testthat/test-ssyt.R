test_that("count_ssytx", {
  N <- count_ssytx(c(4, 3, 3, 2), 5)
  expect_equal(N, 450)
})

test_that("all_ssytx", {
  ssytx <- all_ssytx(c(2, 1), 3)
  expect_length(ssytx, 8L)
})