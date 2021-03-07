test_that("pca_mfd works well", {
  mfdobj <- data_sim_mfd()
  expect_is(pca_mfd(mfdobj), "pca_mfd")
  # Error with one observation
  expect_error(pca_mfd(mfdobj[1]),
               "There is only one observation in the data set")
  # Works with only one variable.
  expect_is(pca_mfd(mfdobj[, 1]), "pca_mfd")
})


