knitr::knit("vignettes/capezza2020.Rmd.orig", "vignettes/capezza2020.Rmd")
knitr::knit("vignettes/centofanti2021.Rmd.orig", "vignettes/centofanti2021.Rmd")
knitr::knit("vignettes/mfd.Rmd.orig", "vignettes/mfd.Rmd")

devtools::build_vignettes()


## After building the website with pkgdown, rename files to get the right image paths
files <- list.files("docs/articles", full.names = TRUE)
files_html <- files[grepl(".html", files)]

for (ff in files_html) {
  text <- readLines(ff)
  new_text <- gsub("../../../R%20packages/funcharts/vignettes/", "", text)
  writeLines(new_text, ff)
}




