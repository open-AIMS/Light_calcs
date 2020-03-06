library(roxygen2)
roxygen2::roxygenise(clean = TRUE)



library(devtools)
library(roxygen2)
library(knitr)
library(R.rsp)
library(digest)

devtools::document()


build()


