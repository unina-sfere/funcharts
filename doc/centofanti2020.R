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
mfdobj <- get_mfd_list(air,
                       grid = 1:24,  
                       n_basis = 15,
                       lambda = 1e-2)
mfdobj_y <- mfdobj[, "NO2"]
mfdobj_x <- mfdobj[, fun_covariates]

## -----------------------------------------------------------------------------
rows1 <- 1:300
rows2 <- 301:355
mfdobj_x1 <- mfdobj_x[rows1]
mfdobj_x2 <- mfdobj_x[rows2]
mfdobj_y1 <- mfdobj_y[rows1]
mfdobj_y2 <- mfdobj_y[rows2]

## -----------------------------------------------------------------------------
mod_fof <- fof_pc(mfdobj_y = mfdobj_y1, mfdobj_x = mfdobj_x1)

## -----------------------------------------------------------------------------
plot_bifd(mod_fof$beta_fd)

## -----------------------------------------------------------------------------
frcc_df <- regr_cc_fof(object = mod_fof,
                       mfdobj_y_new = mfdobj_y2,
                       mfdobj_x_new = mfdobj_x2)
plot_control_charts(frcc_df)

## -----------------------------------------------------------------------------
y_hat <- predict_fof_pc(object = mod_fof,
                        mfdobj_y_new = mfdobj_y2,
                        mfdobj_x_new = mfdobj_x2)

## -----------------------------------------------------------------------------
plot_mon(cclist = frcc_df,
         fd_train = mod_fof$residuals,
         fd_test = y_hat$pred_error[54])

## -----------------------------------------------------------------------------
mfd_list <- get_mfd_list_real_time(data_list = air,
                                   grid = 1:24,
                                   n_basis = 15,
                                   lambda = 1e-2,
                                   n_instants = 7,
                                   start = .5)

mfd_list_x1 <- lapply(mfd_list, function(x) x[rows1, fun_covariates])
mfd_list_x2 <- lapply(mfd_list, function(x) x[rows2, fun_covariates])
mfd_list_y1 <- lapply(mfd_list, function(x) x[rows1, "NO2"])
mfd_list_y2 <- lapply(mfd_list, function(x) x[rows2, "NO2"])

## -----------------------------------------------------------------------------
mod_fof_pc_real_time_list <- fof_pc_real_time(
  mfdobj_y_list = mfd_list_y1,
  mfdobj_x_list = mfd_list_x1)

## -----------------------------------------------------------------------------
cc_list_real_time <- regr_cc_fof_real_time(
  mod_list = mod_fof_pc_real_time_list,
  mfdobj_y_new_list = mfd_list_y2,
  mfdobj_x_new_list = mfd_list_x2
)

## -----------------------------------------------------------------------------
plot_control_charts_real_time(cc_list_real_time, id_num = 54)

