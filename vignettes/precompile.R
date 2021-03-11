knitr::knit("vignettes/capezza2020.Rmd.orig", "vignettes/capezza2020.Rmd")
knitr::knit("vignettes/centofanti2020.Rmd.orig", "vignettes/centofanti2020.Rmd")
knitr::knit("vignettes/mfd.Rmd.orig", "vignettes/mfd.Rmd")

devtools::build_vignettes()
