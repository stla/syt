test_that(".GTpatternToTableau", {
  gt <- lapply(list(
    c(5),
    c(5,4),
    c(3,3,2),
    c(3,3,2,1),
    c(3,3,2,1,1),
    c(3,2,1,1,0,0)
  ), as.integer)
  ssyt <- .GTpatternToTableau(gt)
  expect_identical(
    ssyt,
    list(
      c(1L, 1L, 1L, 5L, 5L),
      c(2L, 2L, 3L, 6L),
      c(3L, 4L),
      4L,
      6L
    )
  )
})

test_that("Number of GT patterns is Kostka number", {
  partitions <- partitions::parts(5)
  Knumbers <- apply(partitions, 2L, function(lambda) {
    apply(partitions, 2L, function(mu) {
      KostkaNumber(lambda, mu)
    })
  })
  counts <- apply(partitions, 2L, function(lambda) {
    apply(partitions, 2L, function(mu) {
      length(GelfandTsetlinPatterns(lambda, mu))
    })
  }) 
  expect_true(all(Knumbers == counts))
})

test_that("Number of skew GT patterns", {
  lambda <- c(4, 2, 1); mu <- c(2, 1); n <- 5
  nskewTableaux <- length(all_ssSkewTableaux(lambda, mu, n))
  compos <- partitions::compositions(sum(lambda)-sum(mu), n)
  skewGTpatterns <- apply(compos, 2L, function(w) {
    skewGelfandTsetlinPatterns(lambda, mu, w)
  }, simplify = FALSE)
  expect_true(
    nskewTableaux == sum(lengths(skewGTpatterns))
  )
})