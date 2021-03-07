data("air")
air <- lapply(air, function(x) x[1:10, , drop = FALSE])
fun_covariates <- names(air)[names(air) != "NO2"]
mfdobj_x <- get_mfd_list(air[fun_covariates], lambda = 1e-2)
y <- rowMeans(air$NO2)

test_that("get_sof_pc_outliers", {
  expect_error(get_sof_pc_outliers(y[1], mfdobj[1:4]))
})
