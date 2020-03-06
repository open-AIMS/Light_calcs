library(devtools)
library(roxygen2)
library(knitr)
library(R.rsp)
library(digest)

roxygen2::roxygenise(clean = TRUE)
build()

# devtools::install_github('AIMS/Light_calcs')
# library(Light_calcs)
