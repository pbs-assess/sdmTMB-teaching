# make clean data frame for sdmTMB course ex

touse <- c("here","tidyverse","sdmTMB") 
lapply(touse, require, character.only=TRUE, quietly=TRUE)

dfc <- readRDS('C://Users//a37789//Documents//coastal-cod//garn-ruse//data//garn_ruse_catch_age.rds')
dfc$year <- as.factor(dfc$year)
dfc$gear2 <- as.factor(dfc$gear2)

# add physical oceanography variables
df1 <- read.table('C://Users//a37789//Documents//coastal-cod//garn-ruse//data//fromJon//Hsig_Hbott_locations_2023.asc', header=TRUE)
dfc <- cbind(dfc, df1[,c("Hsig.ave","Hbott.ave")])
dfc$subarea <- as.factor(dfc$subarea)
dfc$depth2 <- dfc$depth^2

saveRDS(dfc, here('imr-2023','data','garn_ruse_cod_CAA.rds'))

# individual age/length/weight data
# interpolated weight from trawl survey length-weight relationship
# see C:/Users/a37789/Documents/coastal-cod/garn-ruse/code/5-index.R
#   garn_ruse_alk_expanded_LW.rds
