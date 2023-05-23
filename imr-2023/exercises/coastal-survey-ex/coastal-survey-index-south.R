# spatiotemporal index standardization
# coastal trawl survey, swept area biomass of cod
#  note: simplified version of how we use these data

# sdmTMB vignette
# https://pbs-assess.github.io/sdmTMB/index.html

touse <- c("here","tidyverse","sdmTMB","sf","viridis","ggOceanMaps") 
lapply(touse, require, character.only=TRUE, quietly=TRUE)
plotdir <- here('exercises','coastal-survey-ex','plots','south')
resdir <- here('exercises','coastal-survey-ex','results','south')
if(!dir.exists(plotdir)) dir.create(plotdir, recursive=TRUE)
if(!dir.exists(resdir)) dir.create(resdir, recursive=TRUE)

# read clean data frame with cod swept area density (kg / nm2)
# dfc <- readRDS(here('data','survey_catch_clean.rds'))
dfc <- readRDS(here('data','survey_catch_clean_stox.rds'))
dfc <- filter(dfc, !is.na(BottomDepth))
dfc$fyear <- factor(dfc$year)

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

# qqplot from fit5 looks best, spatial + depth

# results from spatial, fit5
# spatial residuals
png(file.path(plotdir, 'resids_spatial.png'), units='in', width=10, height=6, res=150)
print(ggplot(dfc_south, aes(X, Y, col = resids)) + scale_colour_gradient2() +
        geom_point() + facet_wrap(~year) + coord_fixed() + theme_bw())
dev.off()

# prediction grid
res <- 10 # in km
dfc_south$X_bin <- floor(dfc_south$X / res) * res
dfc_south$Y_bin <- floor(dfc_south$Y / res) * res
n_cells <- length(unique(paste0(dfc_south$X_bin, dfc_south$Y_bin)))
dfc_south$XY_bin <- paste(dfc_south$X_bin, dfc_south$Y_bin)
# get mean depth by bin
depth <- dplyr::group_by(dfc_south, XY_bin) %>%
  dplyr::summarise(BottomDepth = mean(BottomDepth))
pred_grid <- data.frame(XY_bin = unique(dfc_south$XY_bin))
pred_grid <- dplyr::left_join(pred_grid, depth)
pred_grid$X <- as.numeric(unlist(lapply(strsplit(pred_grid$XY_bin," "), getElement, 1)))
pred_grid$Y <- as.numeric(unlist(lapply(strsplit(pred_grid$XY_bin," "), getElement, 2)))
pred_grid <- dplyr::select(pred_grid, X, Y, BottomDepth)

# convert dataframe to sf object and set CRS 
pred_grid$X_m <- pred_grid$X*1000
pred_grid$Y_m <- pred_grid$Y*1000
df_sf <- st_as_sf(pred_grid, coords = c("X_m", "Y_m"), crs = get_crs(dfc_south, ll_names=c('lon','lat')))
df_geo <- st_transform(df_sf, crs = 4326) # transform to geographic coordinates (longitude/latitude)
coords <- st_coordinates(df_geo)

# Create a new dataframe with the coordinates
df_coords <- data.frame(
  lon = coords[, 1],
  lat = coords[, 2],
  BottomDepth = pred_grid$BottomDepth
)
# check lat lon conversion
ggplot(df_coords, aes(lon, lat, fill=BottomDepth)) +
   geom_point() + theme_bw() 

# need all covariates (bottom depth, year, X, Y)
predgrid <- crossing(pred_grid, fyear=unique(dfc_south$fyear))
predgrid_geo <- crossing(df_coords, fyear=unique(dfc_south$fyear))
colnames(predgrid) <- c('X','Y','BottomDepth','fyear')
preds <- predict(fit5, newdata = predgrid, return_tmb_object = TRUE, type = "response")
predgrid$est <- preds$data$est
predgrid_geo$est <- predgrid$est

# plot predictions
plot_map <- function(dat, column) {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_tile() +
    coord_fixed() + 
    theme_bw()
}
png(file.path(plotdir, 'preds_spatial.png'), units='in', width=7, height=7, res=150)
print(plot_map(filter(predgrid, fyear == 2022), 'est') + scale_fill_distiller(palette='Oranges', direction=1) + ggtitle("Predicted cod density 2022"))
dev.off()

png(file.path(plotdir, 'preds_spatiotemporal.png'), units='in', width=7, height=7, res=150)
print(plot_map(predgrid, 'est') + facet_wrap(vars(fyear)) + scale_fill_viridis_c() + ggtitle("Predicted cod density"))
dev.off()

png(file.path(plotdir, paste0('spatial_field_pos.png')), units='in', width=7, height=7, res=150)
print(plot_map(preds$data, "omega_s2") + ggtitle("Spatial field positive") + scale_fill_gradient2(trans = 'reverse'))
dev.off()

png(file.path(plotdir, paste0('spatial_field_bin.png')), units='in', width=7, height=7, res=150)
print(plot_map(preds$data, "omega_s1") + ggtitle("Spatial field binomial") + scale_fill_gradient2(trans = 'reverse'))
dev.off()

# better map
plot.yrs <- c(2015, 2022)
toplot <- filter(predgrid_geo, fyear %in% plot.yrs)
maxswept <- max(toplot$est)
g <- basemap(limits = c(5, 12, 61.9, 67.1), land.col = "grey95") +
  scale_y_continuous(breaks = seq(62,67,by=1), name='', expand=c(0,0)) +
  scale_x_continuous(breaks = seq(5,12,by=3), name='', expand=c(0,0)) +
  geom_spatial_point(data = toplot, aes(lon, lat, color = est), size=2) +
  scale_color_distiller(palette='Oranges', direction=1) +
  facet_wrap(vars(fyear)) +
  labs(x = "", y = "", size = "Swept area density cod")+
  theme_bw() +
  theme(axis.text=element_text(size = 10))
ggsave(plot = g, filename = file.path(plotdir,'preds_2years.png'), width = 10, height = 7, units="in", dpi = 300)

# index of abundance
index <- get_index(preds, bias_correct = FALSE, level=0.9)
index$year = as.integer(as.character(index$fyear))
g <- ggplot(index, aes(x=year, y=est)) + 
  ylab('Swept area index (cod biomass)')   +
  geom_ribbon(aes(x=year, ymin=lwr, ymax=upr), color=NA, alpha=.15) +
  scale_y_continuous(expand=c(0.01,0.01), limits = function(x) c(0,range(pretty(x))[2])) +
  geom_line() + 
  geom_point() +
  scale_x_continuous(expand=c(0.02,0.02)) +
  theme_bw() +
  theme(legend.position="top", legend.box.margin = ggplot2::margin(0,0,0,0), legend.margin = ggplot2::margin(0,0,0,0))
png(file.path(plotdir, 'index_0.90ci.png'), units='in', width=10, height=6, res=150)
print(g)
dev.off()

# ----------------------------------------------------------------------------
# compare to stox index
head(index)
# fyear      est      lwr      upr   log_est        se year
# 1  2003 46992.26 31799.45 69443.72 10.757738 0.2374277 2003
# 2  2004 31913.91 23746.03 42891.28 10.370797 0.1797282 2004
# 3  2005 21402.87 14893.71 30756.81  9.971280 0.2204368 2005

index.stox <- read.csv(here('data','coastalcod_totalB2plus_tonn_swept_south.csv'))
head(index.stox)
# Year Ton.mean   Ton.5.   Ton.95.   Ton.sd    Ton.cv
# 1 2003 5998.968 3031.116 10708.695 2302.822 0.3838697
# 2 2004 8403.833 4595.142 12910.226 2656.709 0.3161306
# 3 2005 4429.642 2274.687  6738.668 1335.574 0.3015084

col.labs <- c('year','est','lwr','upr','se')
df1 <- index[,col.labs]
df1$label = 'sdmTMB'
df2 <- index.stox[,c('Year','Ton.mean','Ton.5.','Ton.95.','Ton.cv')]
colnames(df2) = col.labs
df2$label = 'Stox'

df <- rbind(df1, df2)
df <- df %>% group_by(label) %>% mutate(est_stand = est / mean(est),
                                         lo_stand = lwr / mean(est),
                                         hi_stand = upr / mean(est))
years <- range(df$year)
png(file.path(plotdir,'indices_standardized.png'), res=300, units='in', height=8, width=7)
print(ggplot(df, aes(x=year, y=est_stand, ymin=lo_stand, ymax=hi_stand)) +
        geom_hline(yintercept=1, linetype=2) +
        geom_ribbon(aes(fill=label), alpha=0.15) +
        geom_line(aes(color=label), linewidth=1) + 
        facet_wrap(vars(label), ncol=1) +
        ylab('Index (standardized to mean)') +
        xlab('Year') +
        scale_x_continuous(breaks=seq(years[1], years[2], by=2), expand = expansion(mult = c(0.01,0.01))) + 
        scale_y_continuous(limits=c(0, NA), expand = expansion(mult = c(0,0))) +
        # coord_cartesian(ylim=c(0,3)) +
        # scale_y_continuous(limits=c(0, 3), expand = expansion(mult = c(0,0.02))) +
        theme_bw() +
        theme(axis.title = element_text(size=14),
              axis.text.y = element_text(size=10),
              axis.text.x = element_text(size=10),
              legend.title = element_blank(),
              legend.position = 'top'))
dev.off()

# plot CV
png(file.path(plotdir,'indices_cv.png'), res=300, units='in', height=6, width=7)
print(ggplot(df, aes(x=year, y=se)) +
        geom_line(aes(color=label), size=1) + 
        ylab('CV') +
        xlab('') +        
        scale_x_continuous(breaks=seq(years[1], years[2], by=2), expand = expansion(mult = c(0.01,0.01))) + 
        scale_y_continuous(limits=c(0, NA), expand = expansion(mult = c(0,0))) +
        theme_bw() +
        theme(axis.title = element_text(size=14),
              axis.text.y = element_text(size=10),
              axis.text.x = element_text(size=10),
              legend.title = element_blank(),
              legend.position = c(.1,.9)))
dev.off()

# correlation
df <- data.frame(year=index$year,
                 sdmTMB=index$est, 
                 Stox=index.stox$Ton.mean)
png(file.path(plotdir,'indices_correlation.png'), res=300, units='in', height=6, width=7)
ggplot(df, aes(x=sdmTMB, y=Stox)) +
  geom_smooth(method='lm') +
  geom_point() +
  theme_bw()
dev.off()

cor(index$est, index.stox$Ton.mean)
# 0.6220375

