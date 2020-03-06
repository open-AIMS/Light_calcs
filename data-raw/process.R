absorptionDat   <-  read.csv("data-raw/absorption_dat.csv")
modelDat        <-  read.csv("data-raw/model_dat.csv")
solarZenithDat  <-  read.csv("data-raw/solar_zenith_dat.csv")

# Data processing code here...

# This should be the last line.
# Note that names are unquoted.
# I like using overwrite = T so everytime I run the script the 
# updated objects are saved, but the default is overwrite = F
usethis::use_data(absorptionDat, modelDat, solarZenithDat, internal = FALSE, overwrite = TRUE)
