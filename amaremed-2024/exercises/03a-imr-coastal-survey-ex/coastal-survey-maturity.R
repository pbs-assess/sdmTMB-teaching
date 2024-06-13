# spatiotemporal maturity/growth coastal cod
# data: coastal trawl survey
#   note: we have more maturity samples from shallow net survey and reference fleet, could also include

touse <- c("here","tidyverse","sdmTMB","sf","rnaturalearth","viridis","ggOceanMaps") 
lapply(touse, require, character.only=TRUE, quietly=TRUE)
plotdir <- here('exercises','coastal-survey-ex','plots')
resdir <- here('exercises','coastal-survey-ex','results')
if(!dir.exists(plotdir)) dir.create(plotdir, recursive=TRUE)
if(!dir.exists(resdir)) dir.create(resdir, recursive=TRUE)

# read clean data frame with biological data on individual cod (length, weight, age, otolithtype, maturity)
dfa <- readRDS(here('data','survey_age.rds'))

# data cleaning for spatial maturity/growth modeling of coastal cod only
dfa <- dfa %>% filter(startyear > 2002, 
                      gearcondition %in% 1:2,
                      # samplequality %in% 1:3,
                      # samplequality %in% c(1,3),
                      samplequality == 1,
                      gear %in% c(3270, 3271, 3293),
                      commonname == "torsk",
                      otolithtype %in% 1:2)
dfa$lat <- dfa$latitudestart 
nalat <- which(is.na(dfa$lat)) 
dfa$lat[nalat] <-  dfa$latitudeend[nalat] # few latitude NAs filled in
dfa$lon <- dfa$longitudestart 
nalon <- which(is.na(dfa$lon))
dfa$lon[nalon] <-  dfa$longitudeend[nalon] # few longitude NAs filled in
dfa$date <-  as.Date(substr(dfa$stationstopdate, 1, nchar(dfa$stationstopdate)-1), '%Y-%m-%d')
nadate <- which(is.na(dfa$date)) # no NA date
dfa$date[nadate] <-  as.Date(substr(dfa$stationstartdate[nadate], 1, nchar(dfa$stationstartdate[nadate])-1), '%Y-%m-%d')
dfa$year <- dfa$startyear
dfa$quarter <- quarters(dfa$date)

# filter lat north of 62N
dfa <- filter(dfa, lat > 62)

# 2 stocks: 62-67N, north of 67N
dfa$stockarea <- 'north of 67N'
dfa$stockarea[dfa$lat < 67] = '62-67N'

# statistical catch area as factor north -> south
dfa$area <- factor(dfa$area, levels=c(3,4,5,0,6,7))

# depth = average bottom depth
dfa <- dfa %>% rowwise() %>% mutate(depth = mean(c(bottomdepthstart, bottomdepthstop), na.rm=TRUE)) %>% as.data.frame()

# get length in cm and weight in kg to be clear
dfa$length_cm <- 100*dfa$length
dfa$weight_kg <- dfa$individualweight 

# sex: 1 = female, 2 = male, 3 = intersex, some NAs
dfa$sex[dfa$sex == 3] = NA
dfa$sex <- factor(dfa$sex, levels=1:2, labels=c('female','male'))

# maturity data, define maturity as in assessments
#   Q1-2: immature = 1, mature = 2:3, remove 4
#   Q3-4: immature = 1, mature = 2:4
dfa.mat <- dplyr::filter(dfa, !is.na(maturationstage))
torm <- dfa.mat$quarter %in% c('Q1','Q2') & dfa.mat$maturationstage == 4
dfa.mat <- filter(dfa.mat, !torm)
dfa.mat$mat = 0
dfa.mat$mat[dfa.mat$quarter %in% c('Q3','Q4') & dfa.mat$maturationstage %in% 2:4] = 1
dfa.mat$mat[dfa.mat$quarter %in% c('Q1','Q2') & dfa.mat$maturationstage %in% 2:3] = 1 # coastal survey only samples in Q3-Q4, this is for ref fleet data

# growth data
dfa.growth <- dplyr::filter(dfa, !is.na(length_cm) & !is.na(weight_kg)) # all otolith typed have length + weight

# keep columns we'll use
dfa.mat <- dfa.mat %>% select(serialnumber, year, date, quarter, lat, lon, depth, area, 
                      stockarea, mat, maturationstage, age, length_cm, weight_kg, sex)
dfa.growth <- dfa.growth %>% select(serialnumber, year, date, quarter, lat, lon, depth, area, 
                              stockarea, maturationstage, age, length_cm, weight_kg, sex)

dim(dfa.mat)
# [1] 27656    14
dim(dfa.growth)
# [1] 32200    13

saveRDS(dfa.mat, here('data','survey_age_maturity.rds'))
saveRDS(dfa.growth, here('data','survey_age_growth.rds'))

# -----------------------------------------------------------------------------
# map of spatiotemporal sample coverage by year
dfa.mat <- dfa.mat %>% left_join(dfa.mat %>% group_by(year) %>% summarise(N=n())) %>%
  mutate(Label=paste0(year,' (n = ',N,')'))

g <- basemap(limits = c(5, 22, 61.9, 72.5), land.col = "grey95") +
  scale_y_continuous(breaks = seq(62,73,by=1), name='', expand=c(0,0)) +
  scale_x_continuous(breaks = seq(5,22,by=5), name='', expand=c(0,0)) +
  geom_spatial_point(data=dfa.mat, aes(lon, lat, color=stockarea), alpha = 0.6)+
  facet_wrap(vars(Label), nrow=3) +
  theme_bw() +
  theme(axis.text=element_text(size = 10))
ggsave(plot = g, filename = file.path(plotdir,'map_CC_maturity_samples.png'), width = 45, height = 25, units="cm", dpi = 300)

dfa.growth <- dfa.growth %>% left_join(dfa.growth %>% group_by(year) %>% summarise(N=n())) %>%
  mutate(Label=paste0(year,' (n = ',N,')'))

g <- basemap(limits = c(5, 22, 61.9, 72.5), land.col = "grey95") +
  scale_y_continuous(breaks = seq(62,73,by=1), name='', expand=c(0,0)) +
  scale_x_continuous(breaks = seq(5,22,by=5), name='', expand=c(0,0)) +
  geom_spatial_point(data=dfa.growth, aes(lon, lat, color=stockarea), alpha = 0.6)+
  facet_wrap(vars(Label), nrow=3) +
  theme_bw() +
  theme(axis.text=element_text(size = 10))
ggsave(plot = g, filename = file.path(plotdir,'map_CC_growth_samples.png'), width = 45, height = 25, units="cm", dpi = 300)

# -------------------------------------------------------------------------
# big difference in WAA by subarea
WAA2 <- dfa.growth %>% group_by(area, age) %>%
  summarize(meanW_kg = mean(weight_kg),
            se_meanW = sd(weight_kg)/sqrt(n()),
            n = n())
WAA2$lo = WAA2$meanW_kg - 1.96*WAA2$se_meanW
WAA2$hi = WAA2$meanW_kg + 1.96*WAA2$se_meanW
WAA2 <- filter(WAA2, age<13)
ggplot(WAA2, aes(x=age, y=meanW_kg, ymin=lo, ymax=hi, color=area)) +
  geom_pointrange(position=position_dodge(width = .5)) +
  scale_color_brewer(palette = "Dark2", name='') +
  xlab('Age') +
  ylab('Mean weight (kg)') +
  scale_y_continuous(expand=c(0.01,0.01), limits = function(x) c(0,18)) +
  scale_x_continuous(expand=c(0.01,0.01)) +
  theme_bw() +
  theme(strip.text = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  theme(legend.position = c(.07,.75))
ggsave(file.path(plotdir,'waa-survey-area.png'), width = 8, height = 5, units = "in",dpi=200)


