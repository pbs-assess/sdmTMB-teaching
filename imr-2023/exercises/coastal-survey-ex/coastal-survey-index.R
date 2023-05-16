# spatiotemporal index standardization
# coastal trawl survey, swept area biomass of cod
#  note: simplified version of how we use these data

# sdmTMB vignette
# https://pbs-assess.github.io/sdmTMB/index.html

touse <- c("here","tidyverse","sdmTMB","sf","rnaturalearth","viridis","ggOceanMaps") 
lapply(touse, require, character.only=TRUE, quietly=TRUE)
plotdir <- here('exercises','coastal-survey-ex','plots')
resdir <- here('exercises','coastal-survey-ex','results')
if(!dir.exists(plotdir)) dir.create(plotdir, recursive=TRUE)
if(!dir.exists(resdir)) dir.create(resdir, recursive=TRUE)

# read clean data frame with cod swept area density (kg / nm2)
dfc <- readRDS(here('data','survey_catch_clean.rds'))

# -----------------------------------------------------------------------------
# map of spatiotemporal coverage
maxswept <- max(dfc$density)
g <- basemap(limits = c(5, 22, 61.9, 72.5), land.col = "grey95") +
  scale_y_continuous(breaks = seq(62,73,by=1), name='', expand=c(0,0)) +
  scale_x_continuous(breaks = seq(5,22,by=5), name='', expand=c(0,0)) +
  geom_spatial_point(data = dfc %>% filter(density == 0), aes(lon, lat), shape ="x", show.legend=F)+
  geom_spatial_point(data = dfc %>% filter(density > 0), aes(lon, lat, size = density, color=stockarea), alpha = 0.6)+
  facet_wrap(vars(year), nrow=3) +
  scale_size_continuous(limits = c(0,maxswept)) +
  labs(x = "", y = "", size = "Swept area density cod")+
  theme_bw() +
  theme(axis.text=element_text(size = 10))
ggsave(plot = g, filename = file.path(plotdir,'map_swept_density_cod.png'), width = 45, height = 25, units="cm", dpi = 300)



