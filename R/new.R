sim_funcharts <- function() {

  datI <- simulate_mfd(nobs=500)
  datI_tun <- simulate_mfd(nobs=500)
  datII_ic <- simulate_mfd(20)
  datII_oc1 <- simulate_mfd(20,shift_type_x3="A",d_x3=25,d_y_scalar=0.5,shift_type_y="D",d_y=1)
  datII_oc2 <- simulate_mfd(20,shift_type_x3="A",d_x3=50,shift_type_y="D",d_y=2)
  datII_pca <- datII_sof <- datII_fof <- list()
  datII <- list()
  for (ii in 1:4) {
    datII[[ii]] <- rbind(datII_ic[[ii]],datII_oc1[[ii]],datII_oc2[[ii]])
    datII[[ii]] <- rbind(datII_ic[[ii]],datII_oc1[[ii]],datII_oc2[[ii]])
    datII[[ii]] <- rbind(datII_ic[[ii]],datII_oc1[[ii]],datII_oc2[[ii]])
  }
  datII[[5]] <- c(datII_ic[[5]],datII_oc1[[5]],datII_oc2[[5]])
  names(datII) <- names(datI)

  list(
    datI = datI,
    datI_tun = datI_tun,
    datII_pca = datII_pca,
    datII_sof = datII_sof,
    datII_fof = datII_fof
  )
}
