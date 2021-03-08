## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    out.extra='style="display:block; margin: auto"', 
    fig.align="center", 
    tidy=FALSE,
    fig.width=7, 
    fig.height=4)

## -----------------------------------------------------------------------------
library(funcharts)
data("air")
fun_covariates <- names(air)[names(air) != "NO2"]
mfdobj_x <- get_mfd_list(air[fun_covariates],
                       grid = 1:24,  
                       n_basis = 15,
                       lambda = 1e-2)

## -----------------------------------------------------------------------------
y <- rowMeans(air$NO2)

## -----------------------------------------------------------------------------
rows1 <- 1:300
rows2 <- 301:355
mfdobj_x1 <- mfdobj_x[rows1]
mfdobj_x2 <- mfdobj_x[rows2]
y1 <- y[rows1]
y2 <- y[rows2]

## -----------------------------------------------------------------------------
mod <- sof_pc(y = y1, 
              mfdobj_x = mfdobj_x1, 
              selection = "PRESS", 
              single_min_variance_explained = .01)

## -----------------------------------------------------------------------------
plot_mfd(mod$beta)

## -----------------------------------------------------------------------------
plot_bootstrap_sof_pc(mod, nboot = 10)

## -----------------------------------------------------------------------------
cclist_pca <- control_charts_pca(pca = mod$pca, 
                                 components = mod$components, 
                                 newdata = mfdobj_x2)
plot_control_charts(cclist_pca)

## -----------------------------------------------------------------------------
cclist_sof_pc <- control_charts_sof_pc(mod = mod,
                                       y_test = y2,
                                       mfdobj_x_test = mfdobj_x2)
plot_control_charts(cclist_sof_pc)

## -----------------------------------------------------------------------------
ooc_index <- which_ooc(cclist_sof_pc)
ooc_index 

## -----------------------------------------------------------------------------
cont_plot(cclist_sof_pc, 9)

## -----------------------------------------------------------------------------
plot_mon(cclist_sof_pc, mfdobj_x1, mfdobj_x2[9])

## -----------------------------------------------------------------------------
mfd_list <- get_mfd_list_real_time(
  data_list = air[fun_covariates],
  grid = 1:24,
  n_basis = 15,
  lambda = 1e-2,
  n_instants = 7,
  start = .5)
mfd_list1 <- lapply(mfd_list, function(x) x[rows1])
mfd_list2 <- lapply(mfd_list, function(x) x[rows2])

## -----------------------------------------------------------------------------
mod_list <- sof_pc_real_time(y = y1, mfd_real_time_list = mfd_list1)

## -----------------------------------------------------------------------------
cc_list_real_time <- control_charts_sof_pc_real_time(
  mod_list = mod_list,
  y_test = y2,
  mfdobj_x_test = mfd_list2
)

## -----------------------------------------------------------------------------
plot_control_charts_real_time(cc_list_real_time, id_num = 9)

