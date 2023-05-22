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

# --------------------------------------------------------------------------
# south stock
plotdir <- here('exercises','coastal-survey-ex','plots','south')
resdir <- here('exercises','coastal-survey-ex','results','south')
if(!dir.exists(plotdir)) dir.create(plotdir, recursive=TRUE)
if(!dir.exists(resdir)) dir.create(resdir, recursive=TRUE)

dfc_south <- filter(dfc, lat < 67)
dfc_south <- add_utm_columns(dfc_south, c("lon", "lat"), units="km") # UTM zone 32 covers majority of data, 6-12E
mesh <- make_mesh(dfc_south, xy_cols = c("X","Y"), cutoff = 3)
png(file.path(plotdir, 'mesh_south_3km.png'), units='in', width=7, height=7, res=200)
plot(mesh)
dev.off()

# map of spatiotemporal coverage
maxswept <- max(dfc_south$dens_over15cm)
g <- basemap(limits = c(5, 12, 61.9, 67.1), land.col = "grey95") +
  scale_y_continuous(breaks = seq(62,67,by=1), name='', expand=c(0,0)) +
  scale_x_continuous(breaks = seq(5,12,by=3), name='', expand=c(0,0)) +
  geom_spatial_point(data = dfc_south %>% filter(dens_over15cm == 0), aes(lon, lat), shape ="x", show.legend=F)+
  geom_spatial_point(data = dfc_south %>% filter(dens_over15cm > 0), aes(lon, lat, size = dens_over15cm), color='blue', alpha = 0.6)+
  facet_wrap(vars(year), nrow=3) +
  scale_size_continuous(limits = c(0,maxswept)) +
  labs(x = "", y = "", size = "Swept area density cod")+
  theme_bw() +
  theme(axis.text=element_text(size = 10))
ggsave(plot = g, filename = file.path(plotdir,'map_swept_density_cod_stox_south.png'), width = 45, height = 25, units="cm", dpi = 300)


dfc_south$fyear <- factor(dfc_south$year)
fit0 <- sdmTMB(formula = dens_over15cm ~ 0 + fyear,
               data=dfc_south, 
               time="fyear",
               spatial='off',
               spatiotemporal='off',
               mesh=mesh, 
               family = tweedie(link = "log"))

sanity(fit0)

dfc_south$resids <- residuals(fit0)
png(file.path(plotdir, paste0('qqplot_fit0.png')), units='in', width=5, height=5, res=150)
qqnorm(dfc_south$resids)
abline(a = 0, b = 1)
dev.off()

tidy(fit0, "fixed", conf.int=T)
tidy(fit0, "ran_pars", conf.int=T)

# add bottom depth as quadratic
fit1 <- sdmTMB(formula = dens_over15cm ~ 0 + fyear + poly(BottomDepth, 2),
               data=dfc_south, 
               time="fyear",
               spatial='off',
               spatiotemporal='off',
               mesh=mesh, 
               family = tweedie(link = "log"))

sanity(fit1)

dfc_south$resids <- residuals(fit1)
png(file.path(plotdir, paste0('qqplot_fit1.png')), units='in', width=5, height=5, res=150)
qqnorm(dfc_south$resids)
abline(a = 0, b = 1)
dev.off()

tidy(fit1, "fixed", conf.int=T, digits=3) %>% print(n = Inf)
tidy(fit1, "ran_pars", conf.int=T, digits=3) %>% print(n = Inf)

# add spatial
fit2 <- sdmTMB(formula = dens_over15cm ~ 0 + fyear + poly(BottomDepth, 2),
               data=dfc_south, 
               time="fyear",
               spatial='on',
               spatiotemporal='off',
               mesh=mesh, 
               family = tweedie(link = "log"),
               control = sdmTMBcontrol(newton_loops = 1))

sanity(fit2)

dfc_south$resids <- residuals(fit2)
png(file.path(plotdir, paste0('qqplot_fit2.png')), units='in', width=5, height=5, res=150)
qqnorm(dfc_south$resids)
abline(a = 0, b = 1)
dev.off()

tidy(fit2, "fixed", conf.int=T, digits=3) %>% print(n = Inf)
tidy(fit2, "ran_pars", conf.int=T, digits=3) %>% print(n = Inf)

# spatial no depth
fit2b <- sdmTMB(formula = dens_over15cm ~ 0 + fyear,
               data=dfc_south, 
               time="fyear",
               spatial='on',
               spatiotemporal='off',
               mesh=mesh, 
               family = tweedie(link = "log"),
               control = sdmTMBcontrol(newton_loops = 1))

sanity(fit2b)

dfc_south$resids <- residuals(fit2b)
png(file.path(plotdir, paste0('qqplot_fit2b.png')), units='in', width=5, height=5, res=150)
qqnorm(dfc_south$resids)
abline(a = 0, b = 1)
dev.off()

tidy(fit2b, "fixed", conf.int=T, digits=3) %>% print(n = Inf)
tidy(fit2b, "ran_pars", conf.int=T, digits=3) %>% print(n = Inf)

# add spatiotemporal IID
fit3 <- sdmTMB(formula = dens_over15cm ~ 0 + fyear + poly(BottomDepth, 2),
               data=dfc_south, 
               time="fyear",
               spatial='on',
               spatiotemporal='IID',
               mesh=mesh, 
               family = tweedie(link = "log"))
sanity(fit3)

dfc_south$resids <- residuals(fit3)
png(file.path(plotdir, paste0('qqplot_fit3.png')), units='in', width=5, height=5, res=150)
qqnorm(dfc_south$resids)
abline(a = 0, b = 1)
dev.off()

tidy(fit3, "fixed", conf.int=T, digits=3) %>% print(n = Inf)
tidy(fit3, "ran_pars", conf.int=T, digits=3) %>% print(n = Inf)

# add spatiotemporal AR1
# takes awhile
fit4 <- sdmTMB(formula = dens_over15cm ~ 0 + fyear + poly(BottomDepth, 2),
               data=dfc_south, 
               time="year",
               spatial='on',
               spatiotemporal='ar1',
               mesh=mesh, 
               family = tweedie(link = "log"),
               silent=FALSE)
sanity(fit4)

dfc_south$resids <- residuals(fit4)
png(file.path(plotdir, paste0('qqplot_fit4.png')), units='in', width=5, height=5, res=150)
qqnorm(dfc_south$resids)
abline(a = 0, b = 1)
dev.off()

tidy(fit4, "fixed", conf.int=T, digits=3) %>% print(n = Inf)
tidy(fit4, "ran_pars", conf.int=T, digits=3) %>% print(n = Inf)

# add spatial
fit5 <- sdmTMB(formula = dens_over15cm ~ 0 + fyear + poly(BottomDepth, 2),
               data=dfc_south, 
               time="fyear",
               spatial='on',
               spatiotemporal='off',
               mesh=mesh, 
               family = delta_gamma(),
               control = sdmTMBcontrol(newton_loops = 1))
sanity(fit5)

dfc_south$resids <- residuals(fit5)
png(file.path(plotdir, paste0('qqplot_fit5.png')), units='in', width=5, height=5, res=150)
qqnorm(dfc_south$resids)
abline(a = 0, b = 1)
dev.off()

tidy(fit5, "fixed", conf.int=T, digits=3) %>% print(n = Inf)
tidy(fit5, "ran_pars", conf.int=T, digits=3) %>% print(n = Inf)

# add spatial no depth
fit5b <- sdmTMB(formula = dens_over15cm ~ 0 + fyear,
               data=dfc_south, 
               time="fyear",
               spatial='on',
               spatiotemporal='off',
               mesh=mesh, 
               family = delta_gamma(),
               control = sdmTMBcontrol(newton_loops = 1))
sanity(fit5b)

dfc_south$resids <- residuals(fit5b)
png(file.path(plotdir, paste0('qqplot_fit5b.png')), units='in', width=5, height=5, res=150)
qqnorm(dfc_south$resids)
abline(a = 0, b = 1)
dev.off()

tidy(fit5b, "fixed", conf.int=T, digits=3) %>% print(n = Inf)
tidy(fit5b, "ran_pars", conf.int=T, digits=3) %>% print(n = Inf)

# results from spatial, fit5
# spatial residuals
png(file.path(plotdir, 'resids_spatial.png'), units='in', width=10, height=6, res=150)
print(ggplot(dfc_south, aes(X, Y, col = resids)) + scale_colour_gradient2() +
        geom_point() + facet_wrap(~year) + coord_fixed() + theme_bw())
dev.off()

# prediction grid
coordinates(dfc_south) <- ~ X + Y 
res <- raster(dfc_south, resolution=10, vals=TRUE)
r.pts <- as.data.frame(rasterToPoints(res, na.rm=TRUE))

predgrid <- crossing(r.pts, fyear=unique(dfc_south$fyear))
colnames(predgrid) <- c('X','Y','fyear')
preds <- predict(fit5b, newdata = predgrid, return_tmb_object = TRUE, type = "response")
predgrid$est <- preds$data$est

# stepsize <- 5
# nxy <- round(c(diff(range(coords[,1])), diff(range(coords[,2])))/stepsize)
# projgrid = INLA::inla.mesh.projector(mesh$mesh)
# predgrid = as.data.frame(projgrid$lattice$loc)
# names(predgrid) <- c("X","Y")
predgrid <- crossing(predgrid, fyear=unique(dfc_south$fyear))
preds <- predict(fit5b, newdata = predgrid, return_tmb_object = TRUE, type = "response")
predgrid$est <- preds$data$est

index <- get_index(preds, bias_correct = FALSE)

# plot predictions
plot_map <- function(dat, column) {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    coord_fixed() + theme_bw()
}
png(file.path(plotdir, 'preds_spatial.png'), units='in', width=10, height=6, res=150)
print(plot_map(filter(predgrid, fyear == 2022), 'est') + scale_fill_viridis_c() + ggtitle("Predicted cod density"))
dev.off()

png(file.path(plotdir, paste0('spatial_field_pos.png')), units='in', width=10, height=6, res=150)
print(plot_map(preds$data, "omega_s2") + ggtitle("Spatial field positive") + scale_fill_gradient2() + guides(fill = guide_colorbar(reverse = TRUE)))
dev.off()

png(file.path(plotdir, paste0('spatial_field_bin.png')), units='in', width=10, height=6, res=150)
print(plot_map(preds$data, "omega_s1") + ggtitle("Spatial field binomial") + scale_fill_gradient2() + guides(fill = guide_colorbar(reverse = TRUE)))
dev.off()

# index of abundance
index$year = as.integer(as.character(index$fyear))
g <- ggplot(index, aes(x=year, y=est)) + 
  ylab('Swept area index')   +
  geom_ribbon(aes(x=year, ymin=lwr, ymax=upr), color=NA, alpha=.15) +
  scale_y_continuous(expand=c(0.01,0.01), limits = function(x) c(0,range(pretty(x))[2])) +
  geom_line() + 
  geom_point() +
  scale_x_continuous(expand=c(0.02,0.02)) +
  theme_bw() +
  theme(legend.position="top", legend.box.margin = ggplot2::margin(0,0,0,0), legend.margin = ggplot2::margin(0,0,0,0))
png(file.path(plotdir, 'index.png'), units='in', width=10, height=6, res=150)
print(g)
dev.off()
