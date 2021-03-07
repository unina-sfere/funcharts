data("air")
air <- lapply(air, function(x) x[1:10, , drop = FALSE])
fun_covariates <- names(air)[names(air) != "NO2"]
mfdobj <- get_mfd_list(air,
                       grid = 1:24,
                       n_basis = 13,
                       lambda = 1e-2)
mfdobj_y <- mfdobj[, "NO2"]
mfdobj_x <- mfdobj[, fun_covariates]

test_that("fof_pc with one obs", {
  expect_error(fof_pc(mfdobj_y = mfdobj_y[1], mfdobj_x = mfdobj_x[1]),
               "There is only one observation in the data set")
  expect_error(fof_pc(mfdobj_y = mfdobj_y[1:10], mfdobj_x = mfdobj_x[1:2]),
               paste0("mfdobj_y and mfdobj_x must have ",
                      "the same number of observations."))
})

test_that("fof_pc works with multiple responses", {
  expect_is(fof_pc(mfdobj_y = mfdobj_x[1:5], mfdobj_x = mfdobj_x[6:10]),
            "list")
})

test_that("fof_pc", {
  expect_error(
    predict_fof_pc("not_a_list"),
    "object must be a list produced by fof_pc."
  )
  expect_error(
    predict_fof_pc(list("not from sof_pc" = 1)),
    "object must be a list produced by fof_pc."
  )
})
