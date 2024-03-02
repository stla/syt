test_that("Skew tableaux", {
  ssstx <- all_ssSkewTableaux(c(4, 2, 1), c(2, 1), 3)
  expect_length(ssstx, 54L)
  checks <- vapply(ssstx, isSkewTableau, logical(1L))
  expect_true(all(checks))
  dssstx <- lapply(ssstx, dualSkewTableau)
  dchecks <- vapply(dssstx, isSkewTableau, logical(1L))
  expect_true(all(dchecks))
  expect_true(all(vapply(ssstx, isSemistandardSkewTableau, logical(1L))))
})

test_that("Standard skew tableaux", {
  ssstx <- all_ssSkewTableaux(c(4, 2, 1), c(2, 1), 4)
  standard <- Filter(isStandardSkewTableau, ssstx)
  expect_length(standard, 12L)
})