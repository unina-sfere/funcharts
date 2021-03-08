## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    out.extra='style="display:block; margin: auto"', 
    fig.align="center", 
    tidy=FALSE,
    fig.width=7, 
    fig.height=4)

## ---- message = FALSE---------------------------------------------------------
library(funcharts)
data("air")
fun_covariates <- names(air)[names(air) != "NO2"]
mfdobj_x <- get_mfd_list(air[fun_covariates],
                       grid = 1:24,  
                       n_basis = 15,
                       lambda = 1e-2)

## -----------------------------------------------------------------------------
rows1 <- 1:300
rows2 <- 301:355
mfdobj_x1 <- mfdobj_x[rows1]
mfdobj_x2 <- mfdobj_x[rows2]

## -----------------------------------------------------------------------------
y <- rowMeans(air$NO2)
y1 <- y[rows1]
y2 <- y[rows2]

## -----------------------------------------------------------------------------
plot_mfd(mfdobj_x1)
plot_mfd(mfdobj_x1[1:10, c("CO", "C6H6")])

## -----------------------------------------------------------------------------
dat <- data.frame(id = unique(mfdobj_x1$raw$id)) %>% 
  mutate(id_greater_than_100 = as.numeric(id) > 100)
ggplot() + 
  geom_mfd(mapping = aes(col = id_greater_than_100), 
           mfdobj = mfdobj_x1, 
           data = dat, 
           alpha = .2, 
           lwd = .3)

## -----------------------------------------------------------------------------
pca <- pca_mfd(mfdobj_x1)
plot_pca_mfd(pca, harm = 1:3)

## -----------------------------------------------------------------------------
mod <- sof_pc(y = y1, mfdobj_x = mfdobj_x1, selection = "PRESS", single_min_variance_explained = .01)

## -----------------------------------------------------------------------------
plot_bootstrap_sof_pc(mod, nboot = 10, ncores = parallel::detectCores())

## -----------------------------------------------------------------------------
cclist_pca <- control_charts_pca(pca = pca, 
                                 components = mod$components, 
                                 newdata = mfdobj_x2)
plot_control_charts(cclist_pca)

## -----------------------------------------------------------------------------
cclist_sof_pc <- control_charts_sof_pc(
  mod = mod,
  mfdobj_x_test = mfdobj_x2, 
  y_test = y2)
plot_control_charts(cclist_sof_pc)

## -----------------------------------------------------------------------------
ooc_index <- which_ooc(cclist_sof_pc)
ooc_index 

## -----------------------------------------------------------------------------
cont_plot(cclist_sof_pc, 9)
plot_mon(cclist_sof_pc, mfdobj_x1, mfdobj_x2[9])
cont_plot(cclist_sof_pc, 25)
plot_mon(cclist_sof_pc, mfdobj_x1, mfdobj_x2[25])

## -----------------------------------------------------------------------------
mfd_list <- get_mfd_list_real_time(
  data_list = air[fun_covariates],
  grid = 1:24,
  lambda = 1e-2,
  n_instants = 7,
  start = .5
)
mfd_list1 <- lapply(mfd_list, function(x) x[rows1])
mfd_list2 <- lapply(mfd_list, function(x) x[rows2])

mod_list <- sof_pc_real_time(mfd_list1, y1, ncores = parallel::detectCores())

cc_list_real_time <- control_charts_sof_pc_real_time(
  mod_list = mod_list,
  mfdobj_x_test = mfd_list2,
  y_test = y2,
  ncores = parallel::detectCores()
)

plot_control_charts_real_time(cc_list_real_time, id_num = 25)

## -----------------------------------------------------------------------------
library(FRegSigCom)
library(funcharts)
fun_covariates <- names(air)[names(air) != "NO2"]
mfdobj_x <- get_mfd_list(air[fun_covariates],
                         grid = 1:24,
                         n_basis = 15,
                         lambda = 1e-2,
                         ncores = parallel::detectCores())
mfdobj_y <- get_mfd_list(air["NO2"],
                         grid = 1:24,
                         n_basis = 15,
                         lambda = 1e-2,
                         ncores = parallel::detectCores())
rows1 <- 1:300
rows2 <- 301:355
mfdobj_x1 <- mfdobj_x[rows1]
mfdobj_x2 <- mfdobj_x[rows2]
mfdobj_y1 <- mfdobj_y[rows1]
mfdobj_y2 <- mfdobj_y[rows2]

mod_fof <- fof_pc(
  mfdobj_y = mfdobj_y1,
  mfdobj_x = mfdobj_x1)

frcc_df <- regr_cc_fof(
  object = mod_fof,
  mfdobj_y_new = mfdobj_y2,
  mfdobj_x_new = mfdobj_x2)

plot_control_charts(frcc_df)

y_hat <- predict_fof_pc(
  object = mod_fof,
  mfdobj_y_new = mfdobj_y2,
  mfdobj_x_new = mfdobj_x2)

plot_mon(cclist = frcc_df,
         fd_train = mod_fof$residuals,
         fd_test = y_hat$pred_error[28])

mfd_list_x <- get_mfd_list_real_time(
  data_list = air[fun_covariates],
  grid = 1:24,
  lambda = 1e-2,
  n_instants = 7,
  start = .5,
  ncores = parallel::detectCores()
)

mfd_list_y <- get_mfd_list_real_time(
  data_list = air["NO2"],
  grid = 1:24,
  lambda = 1e-2,
  n_instants = 7,
  start = .5,
  ncores = parallel::detectCores()
)

mfd_list_x1 <- lapply(mfd_list_x, function(x) x[rows1])
mfd_list_x2 <- lapply(mfd_list_x, function(x) x[rows2])
mfd_list_y1 <- lapply(mfd_list_y, function(x) x[rows1])
mfd_list_y2 <- lapply(mfd_list_y, function(x) x[rows2])

mod_fof_pc_real_time_list <- fof_pc_real_time(
  mfdobj_y_list = mfd_list_y1,
  mfdobj_x_list = mfd_list_x1,
  ncores = parallel::detectCores())

cc_list_real_time <- regr_cc_fof_real_time(
  mod_list = mod_fof_pc_real_time_list,
  mfdobj_y_new_list = mfd_list_y2,
  mfdobj_x_new_list = mfd_list_x2,
  ncores = parallel::detectCores()
)

plot_control_charts_real_time(cc_list_real_time, id_num = 28)


