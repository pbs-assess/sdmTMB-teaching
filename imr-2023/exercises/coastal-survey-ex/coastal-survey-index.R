# spatiotemporal index standardization
# coastal trawl survey, swept area biomass of cod
#  note: simplified version of how we use these data

# sdmTMB vignette
# https://pbs-assess.github.io/sdmTMB/index.html

touse <- c("here","tidyverse","sdmTMB","sf","rnaturalearth","viridis","ggOceanMaps","alphahull") 
lapply(touse, require, character.only=TRUE, quietly=TRUE)
plotdir <- here('exercises','coastal-survey-ex','plots')
resdir <- here('exercises','coastal-survey-ex','results')
if(!dir.exists(plotdir)) dir.create(plotdir, recursive=TRUE)
if(!dir.exists(resdir)) dir.create(resdir, recursive=TRUE)

# read clean data frame with cod swept area density (kg / nm2)
# dfc <- readRDS(here('data','survey_catch_clean.rds'))
dfc <- readRDS(here('data','survey_catch_clean_stox.rds'))
dfc <- filter(dfc, !is.na(BottomDepth))

# -----------------------------------------------------------------------------
# map of spatiotemporal coverage
maxswept <- max(dfc$dens_over15cm)
g <- basemap(limits = c(5, 22, 61.9, 72.5), land.col = "grey95") +
  scale_y_continuous(breaks = seq(62,73,by=1), name='', expand=c(0,0)) +
  scale_x_continuous(breaks = seq(5,22,by=5), name='', expand=c(0,0)) +
  geom_spatial_point(data = dfc %>% filter(dens_over15cm == 0), aes(lon, lat), shape ="x", show.legend=F)+
  geom_spatial_point(data = dfc %>% filter(dens_over15cm > 0), aes(lon, lat, size = dens_over15cm, color=stockarea), alpha = 0.6)+
  facet_wrap(vars(year), nrow=3) +
  scale_size_continuous(limits = c(0,maxswept)) +
  labs(x = "", y = "", size = "Swept area density cod")+
  theme_bw() +
  theme(axis.text=element_text(size = 10))
ggsave(plot = g, filename = file.path(plotdir,'map_swept_density_cod_stox.png'), width = 45, height = 25, units="cm", dpi = 300)


dfc <- add_utm_columns(dfc, c("lon", "lat"), units="km") # UTM zone 32 covers majority of data, 6-12E
mesh <- make_mesh(dfc, xy_cols = c("X","Y"), cutoff = 3)
png(file.path(plotdir, 'mesh_3km.png'), units='in', width=7, height=7, res=200)
plot(mesh)
dev.off()

dfc$fyear <- factor(dfc$year)
fit0 <- sdmTMB(formula = dens_over15cm ~ 0 + fyear,
               data=dfc, 
               time="fyear",
               spatial='off',
               spatiotemporal='off',
               mesh=mesh, 
               family = tweedie(link = "log"))

sanity(fit0)

dfc$resids <- residuals(fit0)
png(file.path(plotdir, paste0('qqplot_fit0.png')), units='in', width=5, height=5, res=150)
qqnorm(dfc$resids)
abline(a = 0, b = 1)
dev.off()

tidy(fit0, "fixed", conf.int=T)
tidy(fit0, "ran_pars", conf.int=T)

# add bottom depth as quadratic
fit1 <- sdmTMB(formula = dens_over15cm ~ 0 + fyear + poly(BottomDepth, 2),
               data=dfc, 
               time="fyear",
               spatial='off',
               spatiotemporal='off',
               mesh=mesh, 
               family = tweedie(link = "log"))

sanity(fit1)

dfc$resids <- residuals(fit1)
png(file.path(plotdir, paste0('qqplot_fit1.png')), units='in', width=5, height=5, res=150)
qqnorm(dfc$resids)
abline(a = 0, b = 1)
dev.off()

tidy(fit1, "fixed", conf.int=T, digits=3) %>% print(n = Inf)
tidy(fit1, "ran_pars", conf.int=T, digits=3) %>% print(n = Inf)

# add spatial
fit2 <- sdmTMB(formula = dens_over15cm ~ 0 + fyear + poly(BottomDepth, 2),
               data=dfc, 
               time="fyear",
               spatial='on',
               spatiotemporal='off',
               mesh=mesh, 
               family = tweedie(link = "log"))

sanity(fit2)

dfc$resids <- residuals(fit2)
png(file.path(plotdir, paste0('qqplot_fit2.png')), units='in', width=5, height=5, res=150)
qqnorm(dfc$resids)
abline(a = 0, b = 1)
dev.off()

tidy(fit2, "fixed", conf.int=T, digits=3) %>% print(n = Inf)
tidy(fit2, "ran_pars", conf.int=T, digits=3) %>% print(n = Inf)

# add spatiotemporal IID
fit3 <- sdmTMB(formula = dens_over15cm ~ 0 + fyear + poly(BottomDepth, 2),
               data=dfc, 
               time="fyear",
               spatial='on',
               spatiotemporal='IID',
               mesh=mesh, 
               family = tweedie(link = "log"))

sanity(fit3)

dfc$resids <- residuals(fit3)
png(file.path(plotdir, paste0('qqplot_fit3.png')), units='in', width=5, height=5, res=150)
qqnorm(dfc$resids)
abline(a = 0, b = 1)
dev.off()

tidy(fit3, "fixed", conf.int=T, digits=3) %>% print(n = Inf)
tidy(fit3, "ran_pars", conf.int=T, digits=3) %>% print(n = Inf)

