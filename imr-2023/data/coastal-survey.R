# Brian Stock
# May 2023
# clean survey data and calculate cod swept area density (kg / nm2)

# read in catch level data
dfc <- readRDS(here('data','survey_catch.rds'))

# filter / clean data
#   startyear: 2003-, changes in survey design, cod swept area index starts 2003
#   samplequality: 1, 3
#   gearcondition: 1-2
#   gear: 3270, 3271, 3293
#   lat: north of 62N
dfc <- dfc %>% filter(startyear > 2002, 
                      gearcondition %in% 1:2,
                      # samplequality %in% 1:3,
                      # samplequality %in% c(1,3),
                      samplequality == 1,
                      gear %in% c(3270, 3271, 3293))
dfc$lat <- dfc$latitudestart 
nalat <- which(is.na(dfc$lat)) 
dfc$lat[nalat] <-  dfc$latitudeend[nalat] # few latitude NAs filled in
dfc$lon <- dfc$longitudestart 
nalon <- which(is.na(dfc$lon))
dfc$lon[nalon] <-  dfc$longitudeend[nalon] # few longitude NAs filled in
dfc$date <-  as.Date(substr(dfc$stationstopdate, 1, nchar(dfc$stationstopdate)-1), '%Y-%m-%d')
nadate <- which(is.na(dfc$date)) # no NA date
dfc$date[nadate] <-  as.Date(substr(dfc$stationstartdate[nadate], 1, nchar(dfc$stationstartdate[nadate])-1), '%Y-%m-%d')
dfc$year <- dfc$startyear

# filter lat north of 62N
dfc <- filter(dfc, lat > 62)

# 2 stocks: 62-67N, north of 67N
dfc$stockarea <- 'north of 67N'
dfc$stockarea[dfc$lat < 67] = '62-67N'

# depth = average fishing depth, 2 NAs
dfc <- dfc %>% rowwise() %>% mutate(depth = mean(c(bottomdepthstart, bottomdepthstop), na.rm=TRUE)) %>% as.data.frame()

# add haulid
dfc$haulid <- paste(dfc$date, dfc$serialnumber, sep="-")

# swept area
# filter strange distance and door spread values... 
hist(dfc$trawldoorspread, breaks=100)
abline(v=25, col='red', lty=2)
abline(v=100, col='red', lty=2)
dfc <- filter(dfc, trawldoorspread > 25 & trawldoorspread < 100)
hist(dfc$distance, breaks=50)
abline(v=0, col='red', lty=2)
dfc <- filter(dfc, distance > 0)
dfc$sweptarea <- dfc$distance*dfc$trawldoorspread/1852

# add rows with 0 cod catch before filtering on cod
# can also have multiple rows of the same species when different sample types are taken
# need to add these count, weight, and cpue
dfc <- dfc %>% complete(nesting(haulid, year, date, serialnumber, location, station, gear, distance, depth, lat, lon), 
                        commonname, fill=list(catchcount=0, catchweight=0)) %>% as.data.frame()
dfc <- dfc %>% filter(commonname == "torsk") %>% 
  group_by(haulid) %>%
  summarize(catchcount = sum(catchcount, na.rm=TRUE), 
            catchweight = sum(catchweight, na.rm=TRUE), 
            across(-c(catchcount, catchweight))) %>% 
  filter(!duplicated(haulid)) %>% 
  mutate(density = catchweight/sweptarea) %>% as.data.frame()
dfc <- dfc %>% select(haulid, year, date, lat, lon, depth, area, stockarea, commonname,
                      serialnumber, sweptarea, catchcount, catchweight, density)
saveRDS(dfc, here('data','survey_catch_clean.rds'))

