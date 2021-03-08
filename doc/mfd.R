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
                       n_basis = 5,
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

