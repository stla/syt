test_that("Skew tableaux", {
  ssstx <- all_ssSkewTableaux(c(4, 2, 1), c(2, 1), 3)
  expect_length(ssstx, 54L)
})