#' Get out of control observations from control charts
#'
#' @param cclist
#' A \code{data.frame} produced by
#' \code{\link{control_charts_pca}}, \code{\link{control_charts_sof_pc}},
#' \code{\link{regr_cc_fof}}, or \code{\link{regr_cc_sof}}.
#'
#' @return
#' A \code{data.frame} with the same number of rows as cclist,
#' and the same number of columns
#' apart from the columns indicating control chart limits.
#' Each value is TRUE if the corresponding observation is in control
#' and FALSE otherwise.
#'
#' @export
#' @examples
#' library(funcharts)
#' data("air")
#' air <- lapply(air, function(x) x[201:300, , drop = FALSE])
#' fun_covariates <- c("CO", "temperature")
#' mfdobj_x <- get_mfd_list(air[fun_covariates],
#'                          n_basis = 15,
#'                          lambda = 1e-2)
#' y <- rowMeans(air$NO2)
#' y1 <- y[1:60]
#' y2 <- y[91:100]
#' mfdobj_x1 <- mfdobj_x[1:60]
#' mfdobj_x_tuning <- mfdobj_x[61:90]
#' mfdobj_x2 <- mfdobj_x[91:100]
#' mod <- sof_pc(y1, mfdobj_x1)
#' cclist <- control_charts_sof_pc(mod = mod,
#'                                 y_test = y2,
#'                                 mfdobj_x_test = mfdobj_x2,
#'                                 mfdobj_x_tuning = mfdobj_x_tuning)
#' get_ooc(cclist)
#'
#'
get_ooc <- function(cclist) {

  df <- data.frame(
    id = cclist$id,
    T2 = cclist$T2 > cclist$T2_lim,
    spe = cclist$spe > cclist$spe_lim)
  if ("y" %in% names(cclist))
    df$y <-
      cclist$pred_err > cclist$pred_err_sup |
      cclist$pred_err < cclist$pred_err_inf

  df_cont_T2 <- df_cont_spe <- NULL

  if (sum(grepl("contribution_T2", names(cclist))) > 0) {
    cont_T2 <- cclist %>%
      select(contains("contribution_T2")) %>%
      select(- contains("_lim"))
    cont_T2_lim <- cclist %>%
      select(contains("contribution_T2")) %>%
      select(contains("_lim"))
    df_cont_T2 <- as.data.frame(cont_T2 > cont_T2_lim)

  }
  if (sum(grepl("contribution_spe", names(cclist))) > 0) {
    cont_spe <- cclist %>%
      select(contains("contribution_spe")) %>%
      select(- contains("_lim"))
    cont_spe_lim <- cclist %>%
      select(contains("contribution_spe")) %>%
      select(contains("_lim"))
    df_cont_spe <- as.data.frame(cont_spe > cont_spe_lim)
  }

  bind_cols(df, df_cont_T2, df_cont_spe)

}

#' Get the index of the out of control observations from control charts
#'
#' This function returns a list for each control chart and returns
#' the id of all observations that are out of control in that control chart.
#'
#' @param cclist
#' A \code{data.frame} produced by
#' \code{\link{control_charts_sof_pc}}.
#'
#' @return
#' A list of as many \code{data.frame} objects as
#' the control charts in \code{cclist}.
#' Each data frame has two columns, the \code{n} contains
#' an index number giving the observation in the
#' phase II data set, i.e. 1 for the first observation,
#' 2 for the second, and so on,
#' while the \code{id} column contains the id of the observation,
#' which can be general and
#' depends on the specific data set.
#'
#' @export
#' @examples
#' library(funcharts)
#' data("air")
#' air <- lapply(air, function(x) x[201:300, , drop = FALSE])
#' fun_covariates <- c("CO", "temperature")
#' mfdobj_x <- get_mfd_list(air[fun_covariates],
#'                          n_basis = 15,
#'                          lambda = 1e-2)
#' y <- rowMeans(air$NO2)
#' y1 <- y[1:60]
#' y2 <- y[91:100]
#' mfdobj_x1 <- mfdobj_x[1:60]
#' mfdobj_x_tuning <- mfdobj_x[61:90]
#' mfdobj_x2 <- mfdobj_x[91:100]
#' mod <- sof_pc(y1, mfdobj_x1)
#' cclist <- control_charts_sof_pc(mod = mod,
#'                                 y_test = y2,
#'                                 mfdobj_x_test = mfdobj_x2,
#'                                 mfdobj_x_tuning = mfdobj_x_tuning)
#' which_ooc(cclist)
#'
which_ooc <- function(cclist) {

  df <- get_ooc(cclist)

  T2 <- spe <- y <- NULL

  out <- list()

  if (!is.null(df$T2)) {
    out$T2 <- data.frame(n = which(df$T2), id = df$id[df$T2])

  }

  if (!is.null(df$spe)) {
    out$spe <- data.frame(n = which(df$spe), id = df$id[df$spe])
  }

  if (!is.null(df$y)) {
    out$y <- data.frame(n = which(df$y), id = df$id[df$y])
  }

  out

}


#' Produce contribution plots
#'
#' This function produces a contribution plot from
#' functional control charts for
#' a given observation of a phase II data set, using ggplot.
#'
#' @param cclist
#' A \code{data.frame} produced by
#' \code{\link{control_charts_pca}}, \code{\link{control_charts_sof_pc}}
#' \code{\link{regr_cc_fof}}, or \code{\link{regr_cc_sof}}.
#' @param id_num
#' An index number giving the observation in the
#' phase II data set to be plotted, i.e. 1 for the first observation,
#' 2 for the second, and so on.
#' @param which_plot
#' A character vector.
#' Each value indicates which contribution you want to plot:
#'
#' "T2" indicates contribution to the Hotelling's T^2 statistic,
#'
#' "spe" indicates contribution to the squared prediction error statistic.
#'
#' @param print_id
#' A logical value, if TRUE,
#' it prints also the id of the observation
#' in the title of the ggplot.
#' Default is FALSE.
#'
#' @return
#' A ggplot containing the contributions of functional variables to the
#' monitoring statistics.
#' Each plot is a bar plot, with bars corresponding to contribution values and
#' horizontal black segments denoting corresponding (empirical) upper limits.
#' Bars are coloured by red if contributions exceed their limit.
#' @export
#' @examples
#' library(funcharts)
#' data("air")
#' air <- lapply(air, function(x) x[201:300, , drop = FALSE])
#' fun_covariates <- c("CO", "temperature")
#' mfdobj_x <- get_mfd_list(air[fun_covariates],
#'                          n_basis = 15,
#'                          lambda = 1e-2)
#' y <- rowMeans(air$NO2)
#' y1 <- y[1:60]
#' y2 <- y[91:100]
#' mfdobj_x1 <- mfdobj_x[1:60]
#' mfdobj_x_tuning <- mfdobj_x[61:90]
#' mfdobj_x2 <- mfdobj_x[91:100]
#' mod <- sof_pc(y1, mfdobj_x1)
#' cclist <- control_charts_sof_pc(mod = mod,
#'                                 y_test = y2,
#'                                 mfdobj_x_test = mfdobj_x2,
#'                                 mfdobj_x_tuning = mfdobj_x_tuning)
#' get_ooc(cclist)
#' cont_plot(cclist, 3)
#'
cont_plot <- function(cclist,
                      id_num,
                      which_plot = c("T2", "spe"),
                      print_id = FALSE) {

  title <- paste("Observation", id_num)
  if (print_id) title <- paste(title, cclist$id[id_num])

  if (identical(which_plot, "all")) which_plot <- c("T2", "spe")

  fun_covariates <- cclist %>%
    select(contains("contribution_T2")) %>%
    select(- contains("_lim")) %>%
    names() %>%
    gsub("contribution_T2_", "", .)

  df_hot <- df_spe <- NULL

  if ("T2" %in% which_plot) {
    df_hot <- data.frame(
      variables = factor(fun_covariates, levels = fun_covariates),
      contribution = cclist %>%
        select(contains("contribution_T2"),
                      - contains("_lim")) %>%
        slice(id_num) %>%
        as.numeric(),
      limit = cclist %>%
        select(contains("contribution_T2")) %>%
        select(contains("_lim")) %>%
        slice(id_num) %>%
        as.numeric()) %>%
      mutate(ooc = .data$contribution > .data$limit) %>%
      mutate(type = "Contribution to T2") %>%
      mutate(x = 1:length(.data$variables) - 0.35,
             xend = 1:length(.data$variables) + 0.35)
  }

  if ("spe" %in% which_plot) {
    df_spe <- data.frame(
      variables = factor(fun_covariates, levels = fun_covariates),
      contribution = cclist %>%
        select(contains("contribution_spe"),
                      - contains("_lim")) %>%
        slice(id_num) %>%
        as.numeric(),
      limit = cclist %>%
        select(contains("contribution_spe")) %>%
        select(contains("_lim")) %>%
        slice(id_num) %>%
        as.numeric()) %>%
      mutate(ooc = .data$contribution > .data$limit) %>%
      mutate(type = "Contribution to SPE") %>%
      mutate(x = 1:length(.data$variables) - 0.35,
             xend = 1:length(.data$variables) + 0.35)
  }

  df <- bind_rows(df_hot, df_spe) %>%
    mutate(type = factor(.data$type, levels = c("Contribution to T2",
                                                "Contribution to SPE")))
  ggplot(df) +
    geom_col(aes(x = .data$variables,
                 y = .data$contribution,
                 fill = .data$ooc), width = 0.7) +
    theme_bw() + xlab("") + ylab("") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none") +
    scale_fill_manual(values = c("FALSE" = "grey", "TRUE" = "tomato1")) +
    geom_segment(aes(
      x = .data$x,
      xend = .data$xend,
      y = .data$limit,
      yend = .data$limit), col = "black") +
    facet_wrap(~ type, scales = "free_y", ncol = 1,
               strip.position = "left") +
    theme(strip.background = element_blank(),
          strip.placement = "outside")

}

#' Plot multivariate functional object over the training data set
#'
#' This function plots selected functions in
#' a phase_II monitoring data set against the
#' corresponding training data set to be compared.
#'
#' @param cclist
#' A \code{data.frame} produced by
#' \code{\link{control_charts_pca}}, \code{\link{control_charts_sof_pc}}
#' \code{\link{regr_cc_fof}}, or \code{\link{regr_cc_sof}}.
#' @param fd_train
#' An object of class \code{mfd} containing
#' the training data set of the functional variables.
#' They are plotted in gray in the background.
#' @param fd_test
#' An object of class \code{mfd} containing
#' the phase II data set of the functional variables to be monitored.
#' They are coloured in black or red on the foreground.
#' @param print_id
#' A logical value, if TRUE,
#' it prints also the id of the observation
#' in the title of the ggplot.
#'
#' @return
#' A ggplot of the multivariate functional data.
#' In particular, the multivariate functional data given in
#' \code{fd_train} are plotted on
#' the background in gray, while the multivariate functional data given in
#' \code{fd_test} are
#' plotted on the foreground, the colour
#' of each curve is black or red depending on if that curve
#' was signal as anomalous by at least a contribution plot.
#' @export
#'
#' @examples
#' library(funcharts)
#' data("air")
#' air <- lapply(air, function(x) x[201:300, , drop = FALSE])
#' fun_covariates <- c("CO", "temperature")
#' mfdobj_x <- get_mfd_list(air[fun_covariates],
#'                          n_basis = 15,
#'                          lambda = 1e-2)
#' y <- rowMeans(air$NO2)
#' y1 <- y[1:60]
#' y2 <- y[91:100]
#' mfdobj_x1 <- mfdobj_x[1:60]
#' mfdobj_x_tuning <- mfdobj_x[61:90]
#' mfdobj_x2 <- mfdobj_x[91:100]
#' mod <- sof_pc(y1, mfdobj_x1)
#' cclist <- control_charts_sof_pc(mod = mod,
#'                                 y_test = y2,
#'                                 mfdobj_x_test = mfdobj_x2,
#'                                 mfdobj_x_tuning = mfdobj_x_tuning)
#' plot_control_charts(cclist)
#' cont_plot(cclist, 3)
#' plot_mon(cclist, fd_train = mfdobj_x1, fd_test = mfdobj_x2[3])
#'
plot_mon <- function(cclist, fd_train, fd_test, print_id = FALSE) {

  id_num <- if (length(fd_test$fdnames[[2]]) == 1) {
    id_num <- which(cclist$id == fd_test$fdnames[[2]])
    title <- paste("Observation", id_num)
    if (print_id) title <- paste(title, cclist$id[id_num])
  } else {
    title <- ""
  }

  ooc <- get_ooc(cclist)
  mapping <- NULL
  df <- data.frame(id = cclist$id)
  if (sum(grepl("contribution_T2", names(ooc))) > 0 &
      sum(grepl("contribution_spe", names(ooc)))) {
    cont_T2 <- select(
      ooc,
      contains("contribution_T2"),
      - contains("_lim"))
    cont_spe <- select(
      ooc,
      contains("contribution_spe"),
      - contains("_lim"))
    cont_out <- data.frame(cont_T2 | cont_spe)
    names(cont_out) <- gsub("contribution_T2_", "", names(cont_T2))
    df <- bind_cols(df, cont_out) %>%
      pivot_longer(- id, values_to = "contribution ooc", names_to = "var")
  } else {
    df <- df %>%
      mutate(`contribution ooc` =
               cclist$T2 > cclist$T2_lim |
               cclist$spe > cclist$spe_lim)
  }


  ggplot() +
    geom_mfd(mfdobj = fd_train, col = "grey") +
    geom_mfd(mfdobj = fd_test, data = df,
             mapping = aes(col = .data$`contribution ooc`), lwd = 1) +
    scale_color_manual(values = c("TRUE" = "tomato1", "FALSE" = "black")) +
    ggtitle(title) +
    theme(legend.position = "none")

}
