## code to prepare `crs.pk` dataset goes here

setwd("C:/Users/MAITYA02/OneDrive - Pfizer/Documents/2022/Activity/BLRM/BayesOrdCRM")
crs.pk <- read.csv("BCMA_PK.csv", header=T)[, -1]

setwd("C:/Users/MAITYA02/OneDrive - Pfizer/Documents/2022/Activity/BLRM/BayesOrdCRM/BayesOrdCRM")
usethis::use_data(crs.pk, overwrite = TRUE)
usethis::use_data_raw()