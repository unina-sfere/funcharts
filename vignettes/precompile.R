## Create .Rmd from .Rmd.orig
knitr::knit("vignettes/colosimo2010.Rmd.orig", "vignettes/colosimo2010.Rmd")
knitr::knit("vignettes/capezza2020.Rmd.orig", "vignettes/capezza2020.Rmd")
knitr::knit("vignettes/centofanti2021.Rmd.orig", "vignettes/centofanti2021.Rmd")
knitr::knit("vignettes/mfd.Rmd.orig", "vignettes/mfd.Rmd")

## Move png files generated from .Rmd to vignettes/ folder
file_png <- list.files()
file_png <- file_png[grepl(".png", file_png)]
new_file_png <- paste0("vignettes/", file_png)
file.rename(file_png, new_file_png)

## Build vignettes
devtools::build_vignettes()

## Build website
pkgdown::build_site()

## After building the website with pkgdown,
## rename files to get the right image paths
files <- list.files("docs/articles", full.names = TRUE)
files_html <- files[grepl(".html", files)]

for (ff in files_html) {
  text <- readLines(ff)
  new_text <- gsub("../../../R%20packages/funcharts/articles/", "", text)
  writeLines(new_text, ff)
}




