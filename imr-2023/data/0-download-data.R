# Brian Stock
# May 2023

# Get IMR biological survey data 1970-
#  1. **Download** XML files from biotic database (IMR-wide, yearly)
#  2. **Parse** XML into RDS, at catch and age level, still IMR-wide and yearly files.
#  3. **Filter/aggregate** the yearly IMR-wide RDS into dataset/survey-specific RDS files.

# Datasets:
#   Garn ruse (shallow water, fyke/trammel net), 2013-
#   Coastal survey, cruise series 23, Varanger Stad NOR coastal cruise in autumn, 1995- 
#     trawl index starts 2003 when survey was redesigned, 95-03 otolith typing questionable. 85-94 saithe target

# For each survey, save data at 2 levels:
#   catch   haul level, ie total catch by species
#   age     individual level, only those with ages

# Save files:
#   garn_ruse_age.rds
#   garn_ruse_catch.rds 
#   survey_age.rds 
#   survey_catch.rds

# Years are hard-coded because I have not been able to figure out how to filter by cruise series...

# This is what I (Brian) currently do but should be improved in the future... Feedback from Edvin:
# I (Edvin) will include these in a vignette.
#   https://github.com/StoXProject/RstoxData/issues/294
# If this works for now, you may not bother changing the current version. But do check pt 4. below.
# 1. I recommend retiring the biotic.R script and use the function RstoxData::ReadBiotic in stead. It has some slight changes in naming conventions, and data formatting (data.table in stead of data.frame, and an additional list-nesting to support reading several files), but I found it easy to do the transition. It is also much faster than the biotic.R parser.
# 2. I agree that cruise-numbers are best for accessing or filtering cruises. Unfortunately fisheries dependent data are not assigned cruise-numbers. If you prefer to access by cruise-series and year, the necessary information to link can be found in the cruiseseries-table (https://referenceeditor.hi.no/apps/referenceeditor/v2/tables/cruiseseries). You can navigate to each year via the button 'samples', and from there to each leg of the survey with the button 'cruises'. You can get all of these three tables in an excel file if you search for 'cruiseseries' at https://referenceeditor.hi.no/apps/referenceeditor/v2/tables. You can also access the cruise-API, but I have yet to venture there. It is definitely a candidate for the library planned for API access.
# 3. I recommend filtering by missiontype rather than serialnumbers to obtain reference-fleet data. The available mission types are listed in https://referenceeditor.hi.no/apps/referenceeditor/v2/tables/missionType.
# 4. An often over-looked fact about the biotic-data model is that it allows for several age-readings to be recorded for an individual. It is rarely utilised, but it occationaly happens by mistake. If that occurs, the usual merge you are doing will mess up the individual-count. Probably a quick check that the merge does not mess up the row-count suffices. If it doesn't, let me know and I'll explain how to deal with it. When making the vignette, I will also implement a support function to make this easier.
# 
# Other things to consider:
#   - If you find that you want to do filtering before you merge all the tables to one flat table. I recommend the function RstoxData::FilterBiotic. It allows for propagating filters down in the hierarchy, so that you automatically remove individuals from a station, when you remove a station. It also has options for upwards propagation, which is explained in the function documentation. 
#   - in general I recommend consulting the documentation for the biotic data format, which also contains links to relevant reference-lists: https://www.imr.no/formats/nmdbiotic/v3.1/bioticv3_1_en.html

# ------------------------------------------------------------------------------
library(dplyr)
library(here) # relative paths starting at project folder

source(here("data","biotic.R")) # change employee id to yours, eg a37789
xmldir <- file.path("G:/","data","xml") # annual IMR-wide data files, original XML from database
rdsdir <- file.path("G:/","data","rds") # annual IMR-wide data files
outdir <- file.path("G:/","data","rds-series") # aggregated files by survey/dataset

# -------------------------------------------------------------------------------
# 1. Download XML files from biotic database (IMR-wide, yearly)

# # doesn't take too long with ethernet
# years <- 1970:2022
# for(y in 1:length(years)){
#   url <- paste0("http://tomcat7.imr.no/apis/nmdapi/biotic/v3/", years[y], "/cache")
#   destfile <- file.path(xmldir, paste0("alt1_", years[y], ".xml"))
#   f=CFILE(destfile,mode="wb")
#   ret <- download.file(url, destfile=destfile, mode="wb", method="auto")
#   curlPerform(url=url, writedata=f@ref)
#   close(f)
# }

# above didn't work in 2023
# alternatively, download individual XML files
# https://datasetexplorer.hi.no/apps/datasetexplorer/v2/tools/1

# -----------------------------------------------------------------------------
# 2. Parse XML into RDS, at catch and age level, still IMR-wide and yearly files.
# ** this is the time consuming part **

# Edvin: recommend retiring the biotic.R script and use the function RstoxData::ReadBiotic instead.
# It has some slight changes in naming conventions, and data formatting 
# (data.table in stead of data.frame, and an additional list-nesting to support 
# reading several files), but I found it easy to do the transition. It is also 
# much faster than the biotic.R parser.
years <- 1970:2022
for(y in 1:length(years)){
  source(here("data","biotic.R"))
  destfile <- file.path(xmldir, paste0("alt1_", years[y], ".xml"))
  alt1970 <- parse_biotic(destfile)
  a <- merge(alt1970$mission, alt1970$fishstation,all=T)
  c1970 <- merge(a,alt1970$catchsample,all=T)          #c = catchsample
  i1970 <- merge(c1970,alt1970$individual,all=T)         #i = individual
  a1970 <- merge(i1970,alt1970$agedetermination,all=T)   #a = agedetermination
  c1970 <- as.data.frame(lapply(c1970, type.convert)) 
  a1970 <- as.data.frame(lapply(a1970, type.convert))
  saveRDS(c1970, file=file.path(rdsdir, paste0("catch", years[y], ".rds")))
  saveRDS(a1970, file=file.path(rdsdir, paste0("age", years[y], ".rds")))
  
  rm(list=setdiff(ls(), c("years","outdir","xmldir","rdsdir")))
  gc()
}

#------------------------------------------------------------------------------------------
# 3. Filter/aggregate the yearly IMR-wide RDS into dataset/survey-specific RDS files.
# Edvin: better/easier to use cruise number

# Garn ruse, shallow water fyke/trammel net survey
# https://datasetexplorer.hi.no/apps/datasetexplorer/v2/navigation/Cruiseseries/29
yrs.garn <- c(2013, 2015:2022)
n.yrs <- length(yrs.garn)
cruisecodes <- c('2013506','2015506','2016506','2017507','2018506','2019506','2020508','2021507','2022508')
tmp.c1 <- readRDS(file.path(rdsdir, paste0("catch", yrs.garn[1], ".rds")))
tmp.c1 <- tmp.c1[which(tmp.c1$cruise == cruisecodes[1]),]
tmp.a1 <- readRDS(file.path(rdsdir, paste0("age", yrs.garn[1], ".rds")))
tmp.a1 <- tmp.a1[which(tmp.a1$cruise == cruisecodes[1]),]
for(y in 2:n.yrs){
  tmp.c2 <- readRDS(file.path(rdsdir, paste0("catch", yrs.garn[y], ".rds")))
  tmp.c2 <- tmp.c2[which(tmp.c2$cruise == cruisecodes[y]),]
  tmp.c1 <- merge(tmp.c1, tmp.c2, all=T)
  
  tmp.a2 <- readRDS(file.path(rdsdir, paste0("age", yrs.garn[y], ".rds")))
  tmp.a2 <- tmp.a2[which(tmp.a2$cruise == cruisecodes[y]),]
  tmp.a1 <- merge(tmp.a1, tmp.a2, all=T)
}
saveRDS(tmp.c1, file=file.path(outdir, "garn_ruse_catch.rds"))
saveRDS(tmp.a1, file=file.path(outdir, "garn_ruse_age.rds"))
rm(list=setdiff(ls(), c("years","outdir","xmldir","rdsdir")))

#------------------------------------------------------------------------------------------
# Coastal acoustic-trawl survey

# attempt to get list of cruise codes for the coastal survey
# Rstox::getNMDinfo(c("cs", "Varanger Stad NOR coastal cruise in autumn"))

# copy paste cruise codes into csv
# https://datasetexplorer.hi.no/apps/datasetexplorer/v2/navigation/Cruiseseries/23
cruisecodes <- read.csv(here('data','coastal-survey-codes.csv'))[,1]
cruisecodes.yrs <- as.numeric(substr(cruisecodes, 1, 4))

# from rds files with all IMR data, extract coastal cruise
tmp.c <- readRDS(file.path(rdsdir, paste0("catch", cruisecodes.yrs[1], ".rds")))
tmp.c <- tmp.c[which(tmp.c$cruise == cruisecodes[1]),]
tmp.a <- readRDS(file.path(rdsdir, paste0("age", cruisecodes.yrs[1], ".rds")))
tmp.a <- tmp.a[which(tmp.a$cruise == cruisecodes[1]),]
for(y in 2:length(cruisecodes)){ # not years
  tmp.c2 <- readRDS(file.path(rdsdir, paste0("catch", cruisecodes.yrs[y], ".rds")))
  tmp.c2 <- tmp.c2[which(tmp.c2$cruise == cruisecodes[y]),]
  tmp.c <- merge(tmp.c, tmp.c2, all=T)
  
  tmp.a2 <- readRDS(file.path(rdsdir, paste0("age", cruisecodes.yrs[y], ".rds")))
  tmp.a2 <- tmp.a2[which(tmp.a2$cruise == cruisecodes[y]),]
  tmp.a <- merge(tmp.a, tmp.a2, all=T)
}
saveRDS(tmp.c, file=file.path(outdir, "survey_catch.rds"))
saveRDS(tmp.a, file=file.path(outdir, "survey_age.rds"))
rm(list=setdiff(ls(), c("years","outdir","xmldir","rdsdir")))
