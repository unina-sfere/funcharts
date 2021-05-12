data("air")
air <- lapply(air, function(x) x[1:220, , drop = FALSE])
fun_covariates <- c("CO", "temperature")
mfdobj_x <- get_mfd_list(air[fun_covariates],
                         n_basis = 15,
                         lambda = 1e-2)
mfdobj_y <- get_mfd_list(air["NO2"],
                         n_basis = 15,
                         lambda = 1e-2)
y_scalar <- rowMeans(air$NO2)
y1_scalar <- y_scalar[1:200]
y2_scalar <- y_scalar[201:220]
mfdobj_y1 <- mfdobj_y[1:200]
mfdobj_y2 <- mfdobj_y[201:220]
mfdobj_x1 <- mfdobj_x[1:200]
mfdobj_x2 <- mfdobj_x[201:220]
mod_sof <- sof_pc(y1_scalar, mfdobj_x1)
mod_fof <- fof_pc(mfdobj_y1, mfdobj_x1)

test_that("regr_sof_pc works", {
  expect_error(regr_cc_sof(object = 0,
                           y_new = y2_scalar,
                           mfdobj_x_new = mfdobj_x2),
               "object must be a list produced by sof_pc.")
  expect_error(regr_cc_sof(object = list(123),
                           y_new = y2_scalar,
                           mfdobj_x_new = mfdobj_x2),
               "object must be a list produced by sof_pc.")
  expect_is(regr_cc_sof(object = mod_sof,
                           y_new = y2_scalar,
                           mfdobj_x_new = mfdobj_x2),
               "data.frame")
})
