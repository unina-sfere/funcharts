test_that("AMFEWMA contributions reproduce the displayed statistic", {
  set.seed(11)
  basis <- fda::create.bspline.basis(c(0, 1), nbasis = 7)
  coefficients <- array(rnorm(7 * 8 * 2), c(7, 8, 2))
  data <- funcharts::mfd(coef = coefficients, basisobj = basis,
                         fdnames = list("time", as.character(1:8), c("x", "y")))
  grid <- seq(0, 1, length.out = 7)
  vectors <- qr.Q(qr(matrix(rnorm(14 * 4), nrow = 14)))
  make_model <- function(huber) list(mod_1 = list(
    mfdobj_tuning = data, grid_points = grid, mean_mfdobj = fda::mean.fd(data),
    vectors = vectors, values = c(2, 1.5, 1, .5), lambda = .2,
    k = rep(2, 14), h = 5, huber = huber))

  for (huber in c(TRUE, FALSE)) for (burn in c(0L, 5L)) for (nseq in c(1L, 3L)) {
    model <- make_model(huber)
    set.seed(23)
    with_cont <- funcharts::AMFEWMA_PhaseII(data, model, n_seq_2 = nseq,
                                             l_seq_2 = 20, n_skip = burn)
    rng <- .Random.seed
    set.seed(23)
    without_cont <- funcharts::AMFEWMA_PhaseII(data, model, n_seq_2 = nseq,
                          l_seq_2 = 20, n_skip = burn, include_contributions = FALSE)
    expect_identical(.Random.seed, rng)
    expect_identical(with_cont$RL, without_cont$RL)
    expect_identical(with_cont$V2, without_cont$V2)
    expect_identical(with_cont$cc[, 1:3], without_cont$cc)
    expect_equal(rowSums(with_cont$cc[, c("contribution_T2_x", "contribution_T2_y")]),
                 with_cont$cc$amfewma_monitoring_statistic, tolerance = 1e-9)
    graph <- funcharts::cont_plot(with_cont$cc, 1)
    expect_s3_class(graph, "ggplot")
    expect_true(all(is.na(graph$data$limit)))
    expect_true(all(!graph$data$ooc))
  }
})

test_that("contribution plots retain support for existing T2 and SPE limits", {
  cc <- data.frame(id = "a", contribution_T2_x = 2, contribution_T2_y = 1,
                   contribution_spe_x = .1, contribution_spe_y = .5,
                   contribution_T2_x_lim = 1.5, contribution_T2_y_lim = 2,
                   contribution_spe_x_lim = .2, contribution_spe_y_lim = .4)
  graph <- funcharts::cont_plot(cc, 1)
  expect_equal(graph$data$limit, c(1.5, 2, .2, .4))
  expect_identical(graph$data$ooc, c(TRUE, FALSE, FALSE, TRUE))
  expect_s3_class(graph, "ggplot")
})
