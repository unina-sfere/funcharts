x <- seq(1, 10, length = 25)
y11 <- cos(x)
y21 <- cos(2 * x)
y12 <- sin(x)
y22 <- sin(2 * x)
df <- data.frame(id = factor(rep(1:2, each = length(x))),
                 x = rep(x, times = 2),
                 y1 = c(y11, y21),
                 y2 = c(y12, y22))

data_list <- list(y1 = rbind(y11, y21),
                  y2 = rbind(y12, y22))

data_array <- aperm(simplify2array(data_list), c(2, 1, 3))


test_that("domain must be a vector of 2 numbers", {
  expect_error(get_mfd_df(dt = df,
                          domain = c(1, 10, 20),
                          arg = "x",
                          id = "id",
                          variables = c("y1", "y2")),
               "domain must be a vector with two numbers.")
})

test_that("get_mfd_fd correctly converts fd objects", {
  expect_equal(get_mfd_fd(fd()),
               mfd(coef = array(c(0, 0), dim = c(2, 1, 1)),
                   basisobj = fd()$basis,
                   fdnames = fd()$fdnames))
  expect_equal({
    mfdobj <- data_sim_mfd()
    fdobj <- fd(mfdobj$coefs, mfdobj$basis, mfdobj$fdnames)
    get_mfd_fd(fdobj[1, 1:2])
  },
  mfdobj[1, 1:2])
  expect_equal({
    mfdobj <- data_sim_mfd()
    fdobj <- fd(mfdobj$coefs, mfdobj$basis, mfdobj$fdnames)
    get_mfd_fd(fdobj[1:2, 1])
  },
  mfdobj[1:2, 1])
})


test_that("tensor_product_mfd works with multivariate objects", {
  mfdobj1 <- data_sim_mfd(nobs = 1, nvar = 3)
  mfdobj2 <- data_sim_mfd(nobs = 1, nvar = 2)
  expect_is(tensor_product_mfd(mfdobj1), "bifd")
  expect_is(tensor_product_mfd(mfdobj1, mfdobj2), "bifd")
  expect_equal({
    tp <- tensor_product_mfd(mfdobj1, mfdobj2)
    dim(tp$coef)
  }, c(4, 4, 1, 3 * 2))
})

test_that("scale_mfd returns error with one single obs", {
  mfdobj1 <- data_sim_mfd()
  mfdobj2 <- data_sim_mfd(nobs = 1, seed = 123)
  expect_error({
    scale_mfd(mfdobj2)
  },
  "There is only one observation in the data set")
  expect_error({
    pca_mfd(mfdobj2)
  },
  "There is only one observation in the data set")
  expect_s3_class({
    mfdobj1_scaled <- scale_mfd(mfdobj1)
    mfdobj2_scaled <- scale_mfd(mfdobj2,
                                center = attr(mfdobj1_scaled, "scaled:center"),
                                scale = attr(mfdobj1_scaled, "scaled:scale"))
    mfdobj2_scaled
  },
  "mfd")
})

