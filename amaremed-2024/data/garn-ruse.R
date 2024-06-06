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

saveRDS(dfc, here('data','garn_ruse_cod_CAA.rds'))

# ------------------------------------------------------------------------
# get total biomass per net (simpler ex)
# each fish has length, use coastal survey L-W to estimate weight by fish, sum by net
# see C:/Users/a37789/Documents/coastal-cod/garn-ruse/code/5-index.R
#     C:/Users/a37789/Documents/coastal-cod/garn-ruse/code/01-prep-data.R
df.ind <- readRDS(here('data','garn_ruse_alk_expanded_LW.rds'))
df.haul <- readRDS(here('data','garn_ruse_catch_clean.rds'))

df.ind$date <-  as.Date(substr(df.ind$stationstartdate, 1, nchar(df.ind$stationstartdate)-1), '%Y-%m-%d')
df.ind$year <- as.numeric(format(df.ind$date,'%Y'))
df.ind$haulid <- paste(df.ind$date, df.ind$serialnumber, sep="-")
n.hauls <- dim(df.haul)[1]

# # if you try to get biomass per net by summing across individual weights you will miss the nets with 0 catch
# length(unique(df.ind$haulid))
# n.hauls

df.haul$biomass_kg <- 0 # total cod biomass per net
for(i in 1:n.hauls){
  cur_haul <- df.ind[which(df.ind$haulid == df.haul$haulid[i]),] # individuals for current haul
  df.haul$biomass_kg[i] <- sum(cur_haul$individualweight)
}
# sum(df.haul$biomass_kg)
# sum(df.ind$individualweight)

df.haul <- df.haul %>% rowwise() %>% mutate(depth = mean(c(fishingdepthmin, fishingdepthmax), na.rm=TRUE)) %>% as.data.frame()
df.haul$fyear <- as.factor(df.haul$year)
df.haul <- df.haul %>% select(haulid, year, date, lat, lon, depth, subarea, site, biomass_kg, numbercod, gear, gear2)

saveRDS(df.haul, here('data','garn_ruse_cod_biomass.rds'))
