test_that("simulation functions work", {
  dat <- sim_funcharts(nobs1 = 10, nobs_tun = 10, nobs2 = 9)
  expect_is(dat, "list")
  dat <- simulate_mfd(nobs = 10, save_beta = TRUE)
  expect_is(dat, "list")
  dat <- simulate_mfd(nobs = 10,
                      shift_type_x1 = "A",
                      shift_type_x2 = "A",
                      shift_type_x3 = "A",
                      shift_type_y = "A",
                      R2 = 0.86,
                      d_y = 1,
                      d_x1 = 1,
                      d_x2 = 1,
                      d_x3 = 1)
  expect_is(dat, "list")
  dat <- simulate_mfd(nobs = 10,
                      shift_type_x1 = "B",
                      shift_type_x2 = "B",
                      shift_type_x3 = "B",
                      shift_type_y = "B",
                      R2 = 0.86,
                      d_y = 1,
                      d_x1 = 1,
                      d_x2 = 1,
                      d_x3 = 1)
  expect_is(dat, "list")
  dat <- simulate_mfd(nobs = 10,
                      shift_type_x1 = "C",
                      shift_type_x2 = "C",
                      shift_type_x3 = "C",
                      shift_type_y = "C",
                      R2 = 0.86,
                      d_y = 1,
                      d_x1 = 1,
                      d_x2 = 1,
                      d_x3 = 1)
  expect_is(dat, "list")
  dat <- simulate_mfd(nobs = 10,
                      shift_type_x1 = "D",
                      shift_type_x2 = "D",
                      shift_type_x3 = "D",
                      shift_type_y = "D",
                      R2 = 0.74,
                      d_y = 1,
                      d_x1 = 1,
                      d_x2 = 1,
                      d_x3 = 1)
  expect_is(dat, "list")
  dat <- simulate_mfd(nobs = 10, save_beta = TRUE)
  expect_is(dat, "list")

})
