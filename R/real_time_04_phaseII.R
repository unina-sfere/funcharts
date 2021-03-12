#' Real-time unsupervised multivariate functional control charts
#'
#' This function produces a list of data frames,
#' each of them is produced by \code{\link{control_charts_pca}}
#' and is needed to plot control charts for monitoring
#' multivariate functional covariates
#' each evolving up to an intermediate domain point.
#'
#' @param pca_list
#' A list of lists produced by \code{\link{pca_mfd_real_time}},
#' containing a list of multivariate functional principal component analysis
#' models estimated
#' on functional data each evolving up to an intermediate domain point.
#' @param components
#' See \code{\link{pca_mfd}}.
#' @param mfdobj_x_test
#' A list created using
#' \code{\link{get_mfd_df_real_time}} or
#' \code{get_mfd_list_real_time}, denoting a list of functional data objects
#' in the phase II monitoring data set,
#' each evolving up to an intermediate domain point,
#' with observations of the multivariate functional data.
#' The length of this list and \code{pca_list} must be equal,
#' and their elements in the same position in the list
#' must correspond to the same intermediate domain point.
#' @param mfdobj_x_tuning
#' A list created using
#' \code{\link{get_mfd_df_real_time}} or
#' \code{get_mfd_list_real_time}, denoting a list of functional data objects
#' in the tuning data set (used to estimate control chart limits),
#' each evolving up to an intermediate domain point,
#' with observations of the multivariate functional data
#' The length of this list and \code{pca_list} must be equal,
#' and their elements in the same position in the list
#' must correspond to the same intermediate domain point.
#' If NULL, the training data, i.e. the functional data
#' in \code{pca_list},
#' are also used as the tuning data set.
#' Default is NULL.
#' @param alpha
#' See \code{\link{control_charts_pca}}.
#' @param limits
#' See \code{\link{control_charts_pca}}.
#' @param seed
#' See \code{\link{control_charts_pca}}.
#' @param nfold
#' See \code{\link{control_charts_pca}}.
#' @param ncores
#' If you want parallelization, give the number of cores/threads
#' to be used when creating objects separately for different instants.
#'
#' @return
#' A list of \code{data.frame}s each
#' produced by \code{\link{control_charts_pca}},
#' corresponding to a given instant.
#' @export
#'
#' @seealso \code{\link{pca_mfd_real_time}}, \code{\link{control_charts_pca}}
#'
#' @examples
#' library(funcharts)
#' data("air")
#' air1 <- lapply(air, function(x) x[1:8, , drop = FALSE])
#' air2 <- lapply(air, function(x) x[9:10, , drop = FALSE])
#' mfdobj_x1_list <- get_mfd_list_real_time(air1[c("CO", "temperature")],
#'                                          n_basis = 15,
#'                                          lambda = 1e-2,
#'                                          k_seq = c(0.5, 1))
#' mfdobj_x2_list <- get_mfd_list_real_time(air2[c("CO", "temperature")],
#'                                          n_basis = 15,
#'                                          lambda = 1e-2,
#'                                          k_seq = c(0.5, 1))
#' pca_list <- pca_mfd_real_time(mfdobj_x1_list)
#' cclist <- control_charts_pca_mfd_real_time(
#'   pca_list = pca_list,
#'   components = 1:3,
#'   mfdobj_x_test = mfdobj_x2_list)
#' plot_control_charts_real_time(cclist, 1)
#'
control_charts_pca_mfd_real_time <- function(pca_list,
                                             components,
                                             mfdobj_x_test,
                                             mfdobj_x_tuning = NULL,
                                             alpha = list(T2 = .0125,
                                                          spe = .0125),
                                             limits = "standard",
                                             seed = 0,
                                             nfold = NULL,
                                             ncores = 1) {

  if (is.null(mfdobj_x_tuning)) {
    mfdobj_x_tuning <- lapply(pca_list, function(pca_ii) pca_ii$data)
  }

  kk_seq <- as.numeric(names(pca_list))

  single_k <- function(ii) {

    cclist_kk <- control_charts_pca(
      pca = pca_list[[ii]],
      components = components,
      tuning_data = mfdobj_x_tuning[[ii]],
      newdata = mfdobj_x_test[[ii]],
      alpha = alpha,
      limits = limits,
      seed = seed,
      nfold = nfold,
      ncores = 1
    )

    cclist_kk$arg <- pca_list[[ii]]$data$fdnames[[1]]
    cclist_kk$kk <- kk_seq[ii]
    cclist_kk$range1 <- pca_list[[ii]]$data$basis$rangeval[1]
    cclist_kk$range2 <- pca_list[[ii]]$data$basis$rangeval[2]

    rownames(cclist_kk) <- NULL
    cclist_kk

  }
  if (ncores == 1) {
    control_charts_out <- lapply(seq_along(pca_list), single_k)
  } else {
    if (.Platform$OS.type == "unix") {
      control_charts_out <- mclapply(seq_along(pca_list),
                                     single_k,
                                     mc.cores = ncores)
    } else {
      cl <- makeCluster(ncores)
      clusterExport(cl,
                    c("pca_list",
                      "components",
                      "mfdobj_x_test",
                      "mfdobj_x_tuning",
                      "alpha",
                      "limits",
                      "seed",
                      "nfold",
                      "kk_seq"),
                    envir = environment())
      control_charts_out <- parLapply(cl, seq_along(pca_list), single_k)
      stopCluster(cl)
    }
  }
  control_charts_out %>%
    bind_rows()

}



#' Real-time scalar-on-function regression control charts
#'
#' This function produces a list of data frames,
#' each of them is produced by \code{\link{control_charts_sof_pc}}
#' and is needed to plot control charts for monitoring
#' a scalar response variable and
#' multivariate functional covariates
#' each evolving up to an intermediate domain point.
#'
#' @param mod_list
#' A list of lists produced by \code{\link{sof_pc_real_time}},
#' containing a list of scalar-on-function linear regression models estimated
#' on functional data each evolving up to an intermediate domain point.
#' @param y_test
#' A numeric vector containing the observations of
#' the scalar response variable in the phase II monitoring data set.
#' @param mfdobj_x_test
#' A list created using
#' \code{\link{get_mfd_df_real_time}} or
#' \code{get_mfd_list_real_time}, denoting a list of functional data objects
#' in the phase II monitoring data set,
#' each evolving up to an intermediate domain point,
#' with observations of the multivariate functional covariates.
#' The length of this list and \code{mod_list} must be equal,
#' and their elements in the same position in the list
#' must correspond to the same intermediate domain point.
#' @param mfdobj_x_tuning
#' A list created using
#' \code{\link{get_mfd_df_real_time}} or
#' \code{get_mfd_list_real_time}, denoting a list of functional data objects
#' in the tuning data set (used to estimate control chart limits),
#' each evolving up to an intermediate domain point,
#' with observations of the multivariate functional covariates.
#' The length of this list and \code{mod_list} must be equal,
#' and their elements in the same position in the list
#' must correspond to the same intermediate domain point.
#' If NULL, the training data, i.e. the functional covariates
#' in \code{mod_list},
#' are also used as the tuning data set.
#' Default is NULL.
#' @param alpha
#' See \code{\link{control_charts_sof_pc}}.
#' @param limits
#' See \code{\link{control_charts_sof_pc}}.
#' @param seed
#' See \code{\link{control_charts_sof_pc}}.
#' @param nfold
#' See \code{\link{control_charts_sof_pc}}.
#' @param ncores
#' If you want parallelization, give the number of cores/threads
#' to be used when creating objects separately for different instants.
#'
#' @return
#' A list of \code{data.frame}s each
#' produced by \code{\link{control_charts_sof_pc}},
#' corresponding to a given instant.
#' @export
#'
#' @seealso \code{\link{sof_pc_real_time}}, \code{\link{control_charts_sof_pc}}
#'
#' @examples
#' library(funcharts)
#' data("air")
#' air1 <- lapply(air, function(x) x[1:8, , drop = FALSE])
#' air2 <- lapply(air, function(x) x[9:10, , drop = FALSE])
#' mfdobj_x1_list <- get_mfd_list_real_time(air1[c("CO", "temperature")],
#'                                          n_basis = 15,
#'                                          lambda = 1e-2,
#'                                          k_seq = c(0.5, 1))
#' mfdobj_x2_list <- get_mfd_list_real_time(air2[c("CO", "temperature")],
#'                                          n_basis = 15,
#'                                          lambda = 1e-2,
#'                                          k_seq = c(0.5, 1))
#' y1 <- rowMeans(air1$NO2)
#' y2 <- rowMeans(air2$NO2)
#' mod_list <- sof_pc_real_time(y1, mfdobj_x1_list)
#' cclist <- control_charts_sof_pc_real_time(
#'   mod_list = mod_list,
#'   y_test = y2,
#'   mfdobj_x_test = mfdobj_x2_list)
#' plot_control_charts_real_time(cclist, 1)
#'
control_charts_sof_pc_real_time <- function(mod_list,
                                            y_test,
                                            mfdobj_x_test,
                                            mfdobj_x_tuning = NULL,
                                            alpha = list(
                                              T2 = .0125,
                                              spe = .0125,
                                              y = .025),
                                            limits = "standard",
                                            seed = 0,
                                            nfold = NULL,
                                            ncores = 1) {

  if (is.null(mfdobj_x_tuning)) {
    mfdobj_x_tuning <- lapply(mod_list, function(mod_ii) mod_ii$pca$data)
  }

  kk_seq <- as.numeric(names(mod_list))

  single_k <- function(ii) {

    cclist_kk <- control_charts_sof_pc(
      mod = mod_list[[ii]],
      mfdobj_x_test = mfdobj_x_test[[ii]],
      y_test = y_test,
      mfdobj_x_tuning = mfdobj_x_tuning[[ii]],
      alpha = alpha,
      limits = limits,
      seed = seed,
      nfold = nfold,
      ncores = 1
    )

    cclist_kk$arg <- mod_list[[ii]]$pca$data$fdnames[[1]]
    cclist_kk$kk <- kk_seq[ii]
    cclist_kk$range1 <- mod_list[[ii]]$pca$data$basis$rangeval[1]
    cclist_kk$range2 <- mod_list[[ii]]$pca$data$basis$rangeval[2]

    rownames(cclist_kk) <- NULL
    cclist_kk

  }
  if (ncores == 1) {
    control_charts_out <- lapply(seq_along(mod_list), single_k)
  } else {
    if (.Platform$OS.type == "unix") {
      control_charts_out <- mclapply(seq_along(mod_list),
                                     single_k,
                                     mc.cores = ncores)
    } else {
      cl <- makeCluster(ncores)
      clusterExport(cl,
                    c("mod_list",
                      "mfdobj_x_test",
                      "y_test",
                      "mfdobj_x_tuning",
                      "alpha",
                      "limits",
                      "seed",
                      "nfold",
                      "kk_seq"),
                    envir = environment())
      control_charts_out <- parLapply(cl, seq_along(mod_list), single_k)
      stopCluster(cl)
    }
  }
  control_charts_out %>%
    bind_rows()

}



#' Real-time functional regression control chart
#'
#' This function produces a list of data frames,
#' each of them is produced by \code{\link{regr_cc_fof}}
#' and is needed to plot control charts for monitoring
#' a functional response variable and
#' multivariate functional covariates
#' each evolving up to an intermediate domain point.
#'
#' @param mod_list
#' A list of lists produced by \code{\link{fof_pc_real_time}},
#' containing a list of function-on-function linear regression models estimated
#' on functional data each evolving up to an intermediate domain point.
#' @param mfdobj_y_new_list
#' A list created using
#' \code{\link{get_mfd_df_real_time}} or
#' \code{get_mfd_list_real_time}, denoting a list of functional data objects
#' in the phase II monitoring data set,
#' each evolving up to an intermediate domain point,
#' with observations of the functional response variable
#' The length of this list and \code{mod_list} must be equal,
#' and their elements in the same position in the list
#' must correspond to the same intermediate domain point.
#' @param mfdobj_x_new_list
#' A list created using
#' \code{\link{get_mfd_df_real_time}} or
#' \code{get_mfd_list_real_time}, denoting a list of functional data objects
#' in the phase II monitoring data set,
#' each evolving up to an intermediate domain point,
#' with observations of the multivariate functional covariates.
#' The length of this list and \code{mod_list} must be equal,
#' and their elements in the same position in the list
#' must correspond to the same intermediate domain point.
#' @param mfdobj_y_tuning_list
#' A list created using
#' \code{\link{get_mfd_df_real_time}} or
#' \code{get_mfd_list_real_time}, denoting a list of functional data objects
#' in the tuning data set (used to estimate control chart limits),
#' each evolving up to an intermediate domain point,
#' with observations of the functional response variable.
#' The length of this list and \code{mod_list} must be equal,
#' and their elements in the same position in the list
#' must correspond to the same intermediate domain point.
#' If NULL, the training data, i.e. the functional response
#' in \code{mod_list},
#' is also used as the tuning data set.
#' Default is NULL.
#' @param mfdobj_x_tuning_list
#' A list created using
#' \code{\link{get_mfd_df_real_time}} or
#' \code{get_mfd_list_real_time}, denoting a list of functional data objects
#' in the tuning data set (used to estimate control chart limits),
#' each evolving up to an intermediate domain point,
#' with observations of the multivariate functional covariates.
#' The length of this list and \code{mod_list} must be equal,
#' and their elements in the same position in the list
#' must correspond to the same intermediate domain point.
#' If NULL, the training data, i.e. the functional covariates
#' in \code{mod_list},
#' are also used as the tuning data set.
#' Default is NULL.
#' @param alpha
#' See \code{\link{regr_cc_fof}}.
#' @param ncores
#' If you want parallelization, give the number of cores/threads
#' to be used when creating objects separately for different instants.
#'
#' @return
#' A list of \code{data.frame}s each
#' produced by \code{\link{regr_cc_fof}},
#' corresponding to a given instant.
#' @export
#'
#' @seealso \code{\link{fof_pc_real_time}}, \code{\link{regr_cc_fof}}
#' @examples
#' library(funcharts)
#' data("air")
#' air1 <- lapply(air, function(x) x[1:8, , drop = FALSE])
#' air2 <- lapply(air, function(x) x[9:10, , drop = FALSE])
#' mfdobj_x1_list <- get_mfd_list_real_time(air1[c("CO", "temperature")],
#'                                          n_basis = 15,
#'                                          lambda = 1e-2,
#'                                          k_seq = c(0.5, 1))
#' mfdobj_x2_list <- get_mfd_list_real_time(air2[c("CO", "temperature")],
#'                                          n_basis = 15,
#'                                          lambda = 1e-2,
#'                                          k_seq = c(0.5, 1))
#' mfdobj_y1_list <- get_mfd_list_real_time(air1["NO2"],
#'                                          n_basis = 15,
#'                                          lambda = 1e-2,
#'                                          k_seq = c(0.5, 1))
#' mfdobj_y2_list <- get_mfd_list_real_time(air2["NO2"],
#'                                          n_basis = 15,
#'                                          lambda = 1e-2,
#'                                          k_seq = c(0.5, 1))
#' mod_list <- fof_pc_real_time(mfdobj_y1_list, mfdobj_x1_list)
#' cclist <- regr_cc_fof_real_time(
#'   mod_list = mod_list,
#'   mfdobj_y_new_list = mfdobj_y2_list,
#'   mfdobj_x_new_list = mfdobj_x2_list)
#' plot_control_charts_real_time(cclist, 1)
#'
regr_cc_fof_real_time <- function(mod_list,
                                  mfdobj_y_new_list,
                                  mfdobj_x_new_list,
                                  mfdobj_y_tuning_list = NULL,
                                  mfdobj_x_tuning_list = NULL,
                                  alpha = list(T2 = 0.025, spe = 0.025),
                                  ncores = 1) {


  kk_seq <- as.numeric(names(mod_list))

  single_k <- function(ii) {

    if (is.null(mfdobj_y_tuning_list) | is.null(mfdobj_x_tuning_list)) {
      mfdobj_y_tuning_list_ii <- mod_list[[ii]]$pca_y$data
      mfdobj_x_tuning_list_ii <- mod_list[[ii]]$pca_x$data
    }

    regr_cc_fof_kk <- regr_cc_fof(
      object = mod_list[[ii]],
      mfdobj_y_new = mfdobj_y_new_list[[ii]],
      mfdobj_x_new = mfdobj_x_new_list[[ii]],
      mfdobj_y_tuning = mfdobj_y_tuning_list_ii,
      mfdobj_x_tuning = mfdobj_x_tuning_list_ii,
      alpha = alpha
    )

    regr_cc_fof_kk$arg <- mod_list[[ii]]$pca_x$data$fdnames[[1]]
    regr_cc_fof_kk$kk <- kk_seq[ii]
    regr_cc_fof_kk$range1 <- mod_list[[ii]]$pca_x$data$basis$rangeval[1]
    regr_cc_fof_kk$range2 <- mod_list[[ii]]$pca_x$data$basis$rangeval[2]

    regr_cc_fof_kk

  }

  if (ncores == 1) {
    control_charts_out <- lapply(seq_along(mod_list), single_k)
  } else {
    if (.Platform$OS.type == "unix") {
      control_charts_out <- mclapply(seq_along(mod_list),
                                     single_k,
                                     mc.cores = ncores)
    } else {
      cl <- makeCluster(ncores)
      clusterExport(cl,
                    c("mod_list",
                      "mfdobj_y_new_list",
                      "mfdobj_x_new_list",
                      "mfdobj_y_tuning_list",
                      "mfdobj_x_tuning_list",
                      "alpha",
                      "kk_seq"),
                    envir = environment())
      control_charts_out <- parLapply(cl, seq_along(mod_list), single_k)
      stopCluster(cl)
    }
  }
  control_charts_out %>%
    bind_rows()

}


#' Plot real-time control charts
#'
#' This function produces a ggplot
#' with the desired real-time control charts.
#' It takes as input a list of data frames, produced
#' with functions such as
#' \code{\link{regr_cc_fof_real_time}} and
#' \code{\link{control_charts_sof_pc_real_time}},
#' and the id of the observations for which real-time control charts
#' are desired to be plotted.
#' For each control chart, the solid line corresponds to the
#' profile of the monitoring statistic and it is compared against
#' control limits plotted as dashed lines.
#' If a line is outside its limits it is coloured in red.
#'
#' @param cclist
#' A list of data frames, produced
#' with functions such as
#' \code{\link{regr_cc_fof_real_time}} and
#' \code{\link{control_charts_sof_pc_real_time}},
#' @param id_num
#' An index number giving the observation in the
#' phase II data set to be plotted, i.e. 1 for the first observation,
#' 2 for the second, and so on.
#'
#' @return A ggplot with the real-time functional control charts.
#'
#' @details
#' If the line, representing the profile of the
#' monitoring statistic over the functional domain, is out-of-control,
#' then it is coloured in red.
#'
#' @export
#' @seealso \code{\link{regr_cc_fof_real_time}},
#' \code{\link{control_charts_sof_pc_real_time}}
#' @examples
#' library(funcharts)
#' data("air")
#' air1 <- lapply(air, function(x) x[1:8, , drop = FALSE])
#' air2 <- lapply(air, function(x) x[9:10, , drop = FALSE])
#' mfdobj_x1_list <- get_mfd_list_real_time(air1[c("CO", "temperature")],
#'                                          n_basis = 15,
#'                                          lambda = 1e-2,
#'                                          k_seq = c(0.5, 1))
#' mfdobj_x2_list <- get_mfd_list_real_time(air2[c("CO", "temperature")],
#'                                          n_basis = 15,
#'                                          lambda = 1e-2,
#'                                          k_seq = c(0.5, 1))
#' y1 <- rowMeans(air1$NO2)
#' y2 <- rowMeans(air2$NO2)
#' mod_list <- sof_pc_real_time(y1, mfdobj_x1_list)
#' cclist <- control_charts_sof_pc_real_time(
#'   mod_list = mod_list,
#'   y_test = y2,
#'   mfdobj_x_test = mfdobj_x2_list)
#' plot_control_charts_real_time(cclist, 1)
plot_control_charts_real_time <- function(cclist, id_num) {

  xmin <- cclist$range1[1]
  xmax <- max(cclist$range2)

  cclist_id <- filter(cclist, id == unique(cclist$id)[id_num])

  kk_range <- range(cclist_id$kk)
  kk_seq <- seq(kk_range[1], kk_range[2], l = 1000)

  x_axis_tick <- sort(unique(cclist$kk))
  ndigits <- round(max(log10(x_axis_tick))) + 1
  x_axis_tick <- round(x_axis_tick, 3 - ndigits)

  plot_list <- list()

  if ("T2" %in% names(cclist_id)) {
    df_hot <- cclist_id %>%
      select(.data$id, .data$kk, statistic = .data$T2, UCL = .data$T2_lim)
    f_hot <- approxfun(df_hot$kk, df_hot$statistic)
    f_hot_UCL <- approxfun(df_hot$kk, df_hot$UCL)
    df_hot_plot <- data.frame(
      id = df_hot$id[1],
      kk = kk_seq,
      statistic = f_hot(kk_seq),
      UCL = f_hot_UCL(kk_seq)
    ) %>%
      mutate(ooc = .data$statistic > .data$UCL) %>%
      mutate(group = .data$ooc %>%
               rle() %>%
               unclass() %>%
               data.frame() %>%
               mutate(values = 1:n()) %>%
               inverse.rle()) %>%
      ungroup()

    plot_list$p_hot <- df_hot_plot %>%
      rename(T2 = .data$statistic) %>%
      pivot_longer(cols = c(.data$T2, .data$UCL),
                   names_to = "line_type") %>%
      mutate(
        line_type = factor(
          x = .data$line_type,
          levels = c("UCL", "T2", "LCL")),
        group = paste0(.data$group, .data$line_type),
        ooc = case_when(
          line_type == "T2" ~ ooc,
          line_type %in% c("UCL", "LCL") ~ FALSE
        )) %>%
      ggplot() +
      geom_line(aes(x     = .data$kk,
                    y     = .data$value,
                    lty   = .data$line_type,
                    col   = .data$ooc,
                    group = .data$group)) +
      scale_linetype_manual(values = c("T2" = 1, "UCL" = 2, "LCL" = 2)) +
      scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "tomato1")) +
      guides(colour = FALSE) +
      geom_blank(aes(y = 0)) +
      scale_x_continuous(limits = c(xmin, xmax),
                         breaks = x_axis_tick,
                         expand = c(0, 0)) +
      theme_bw() +
      theme(strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      xlab(cclist_id$arg[1]) +
      ylab("") +
      ggtitle(expression(HOTELLING~T^2~CONTROL~CHART)) +
      theme(legend.title = element_blank())
  }

  if ("spe" %in% names(cclist_id)) {
    df_spe <- cclist_id %>%
      select(.data$id, .data$kk, statistic = .data$spe, UCL = .data$spe_lim)
    f_spe <- approxfun(df_spe$kk, df_spe$statistic)
    f_spe_UCL <- approxfun(df_spe$kk, df_spe$UCL)
    df_spe_plot <- data.frame(
      id = df_spe$id[1],
      kk = kk_seq,
      statistic = f_spe(kk_seq),
      UCL = f_spe_UCL(kk_seq)
    ) %>%
      mutate(ooc = .data$statistic > .data$UCL) %>%
      mutate(group = .data$ooc %>%
               rle() %>%
               unclass() %>%
               data.frame() %>%
               mutate(values = 1:n()) %>%
               inverse.rle()) %>%
      ungroup()

    plot_list$p_spe <- df_spe_plot %>%
      rename(SPE = .data$statistic) %>%
      pivot_longer(cols = c(.data$SPE, .data$UCL),
                   names_to = "line_type") %>%
      mutate(
        line_type = factor(
          x = .data$line_type,
          levels = c("UCL", "SPE", "LCL")),
        group = paste0(.data$group, .data$line_type),
        ooc = case_when(
          line_type == "SPE" ~ ooc,
          line_type %in% c("UCL", "LCL") ~ FALSE
        )) %>%
      ggplot() +
      geom_line(aes(x     = .data$kk,
                    y     = .data$value,
                    lty   = .data$line_type,
                    col   = .data$ooc,
                    group = .data$group)) +
      scale_linetype_manual(values = c("SPE" = 1, "UCL" = 2, "LCL" = 2)) +
      scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "tomato1")) +
      guides(colour = FALSE) +
      geom_blank(aes(y = 0)) +
      scale_x_continuous(limits = c(xmin, xmax),
                         breaks = x_axis_tick,
                         expand = c(0, 0)) +
      theme_bw() +
      theme(strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      xlab(cclist_id$arg[1]) +
      ylab("") +
      ggtitle(expression(SPE~CONTROL~CHART)) +
      theme(legend.title = element_blank())
  }

  if ("pred_err" %in% names(cclist_id)) {
    df_y <- cclist_id %>%
      select(.data$id,
             .data$kk,
             statistic = .data$pred_err,
             UCL = .data$pred_err_sup,
             LCL = .data$pred_err_inf)
    if (sum(!is.na(df_y$statistic)) > 0) {
      f_y <- approxfun(df_y$kk, df_y$statistic)
    } else {
      f_y <- function(x) NA
    }
    f_y_UCL <- approxfun(df_y$kk, df_y$UCL)
    f_y_LCL <- approxfun(df_y$kk, df_y$LCL)
    df_y_plot <- data.frame(
      id = df_y$id[1],
      kk = kk_seq,
      statistic = f_y(kk_seq),
      UCL = f_y_UCL(kk_seq),
      LCL = f_y_LCL(kk_seq),
      type = "`REAL-TIME`~REGRESSION~CONTROL~CHART"
    )
    if (sum(!is.na(df_y$statistic)) > 0) {
      df_y_plot <- df_y_plot %>%
        mutate(ooc = .data$statistic > .data$UCL |
                 .data$statistic < .data$LCL) %>%
        mutate(group = .data$ooc %>%
                 rle() %>%
                 unclass() %>%
                 data.frame() %>%
                 mutate(values = 1:n()) %>%
                 inverse.rle()) %>%
        ungroup()
    } else {
      df_y_plot <- df_y_plot %>%
        mutate(ooc = NA) %>%
        mutate(group = 1) %>%
        ungroup()
    }

    plot_list$p_y <- df_y_plot %>%
      rename(`prediction error [t]` = .data$statistic) %>%
      pivot_longer(cols = c(.data$`prediction error [t]`, .data$UCL, .data$LCL),
                   names_to = "line_type") %>%
      mutate(
        line_type = factor(
          x = .data$line_type,
          levels = c("UCL", "prediction error [t]", "LCL")),
        group = paste0(.data$group, .data$line_type),
        ooc = case_when(
          line_type == "prediction error [t]" ~ ooc,
          line_type %in% c("UCL", "LCL") ~ FALSE
        )) %>%
      ggplot() +
      geom_line(aes(x     = .data$kk,
                    y     = .data$value,
                    lty   = .data$line_type,
                    col   = .data$ooc,
                    group = .data$group)) +
      scale_linetype_manual(values = c("prediction error [t]" = 1, "UCL" = 2, "LCL" = 2)) +
      scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "tomato1")) +
      guides(colour = FALSE) +
      geom_blank(aes(y = 0)) +
      scale_x_continuous(limits = c(xmin, xmax),
                         breaks = x_axis_tick,
                         expand = c(0, 0)) +
      theme_bw() +
      theme(strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      xlab(cclist_id$arg[1]) +
      ylab("") +
      ggtitle(expression(REGRESSION~CONTROL~CHART)) +
      theme(legend.title = element_blank())
  }

  p <- patchwork::wrap_plots(plot_list, ncol = 1) &
    theme(legend.justification = 0)
  suppressWarnings(print(p))

}

