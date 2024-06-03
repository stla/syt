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

test_that("Comparison all_ssSkewTableaux with skewTableauxWithGivenShapeAndWeight", {
  lambda <- c(4, 2, 1); mu <- c(2, 1); n <- 5
  ssstx <- all_ssSkewTableaux(lambda, mu, n)
  ssstxMatrix <- orderedMatrix(do.call(
    rbind,
    lapply(ssstx, function(ssst) {
      Filter(Negate(is.na), do.call(c, ssst))
    })
  ))
  compos <- partitions::compositions(sum(lambda)-sum(mu), n)
  ssstx2 <- do.call(
    c, 
    apply(compos, 2L, function(compo) {
      skewTableauxWithGivenShapeAndWeight(lambda, mu, compo)
    }, simplify = FALSE)
  )
  ssstxMatrix2 <- orderedMatrix(do.call(
    rbind,
    lapply(ssstx2, function(ssst) {
      Filter(Negate(is.na), do.call(c, ssst))
    })
  ))
  expect_true(
    all(ssstxMatrix == ssstxMatrix2)
  )
})