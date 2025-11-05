library(R6)


test_that("Kmeans class can be initialized", {
  kmeans_obj <- Kmeans$new()
  expect_true(is.R6(kmeans_obj))
})

test_that("fit method works", {
  kmeans_obj <- Kmeans$new(n_cluster=5)
  expect_true(kmeans_obj$fit() == 5)
})
