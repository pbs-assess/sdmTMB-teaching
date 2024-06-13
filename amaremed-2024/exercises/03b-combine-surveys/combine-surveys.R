# can we combine the two surveys to produce an aggregate index?
#   coastal trawl survey
#   shallow net / garn ruse survey

touse <- c("here","tidyverse","sdmTMB","sf","viridis","ggOceanMaps") 
lapply(touse, require, character.only=TRUE, quietly=TRUE)
plotdir <- here('exercises','combine-surveys','plots')
resdir <- here('exercises','combine-surveys','results')
if(!dir.exists(plotdir)) dir.create(plotdir, recursive=TRUE)
if(!dir.exists(resdir)) dir.create(resdir, recursive=TRUE)

# read trawl survey data
dfc <- readRDS(here('data','survey_catch_clean_stox.rds'))
dfc <- filter(dfc, !is.na(BottomDepth) & lat < 67)
dfc$depth <- dfc$BottomDepth
dfc$dens <- dfc$dens_over15cm
dfc$gear <- 'Trawl'
dfc$survey <- 'Coastal trawl'
dfc <- dfc %>% select(survey, haulid, year, lon, lat, gear, depth, dens)

# read shallow net survey data
dfc2 <- readRDS(here('data','garn_ruse_cod_biomass.rds'))
dfc2 <- filter(dfc2, lat < 67)
dfc2$dens <- dfc2$biomass_kg
dfc2$survey <- 'Shallow net'
dfc2 <- dfc2 %>% select(survey, haulid, year, lon, lat, gear, depth, dens)

df <- rbind(dfc, dfc2)
df <- add_utm_columns(df, c("lon", "lat"), units="km")
df$fyear <- factor(df$year)
df$gear <- factor(df$gear)
df$survey <- factor(df$survey)
mesh <- make_mesh(df, xy_cols = c("X","Y"), cutoff = 3)
png(file.path(plotdir, 'mesh_3km.png'), units='in', width=7, height=7, res=200)
plot(mesh)
dev.off()

# map of spatiotemporal coverage
g <- basemap(limits = c(5, 12, 61.9, 67.1), land.col = "grey95") +
  scale_y_continuous(breaks = seq(62,67,by=1), name='', expand=c(0,0)) +
  scale_x_continuous(breaks = seq(5,12,by=3), name='', expand=c(0,0)) +
  geom_spatial_point(data = df, aes(lon, lat, color=survey), alpha = 0.6)+
  facet_wrap(vars(year), nrow=3) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.text=element_text(size = 10), legend.position='top')
ggsave(plot = g, filename = file.path(plotdir,'map_coverage.png'), width = 14, height = 10, units="in", dpi = 300)

# trawl survey only has one gear, so careful of depth/gear/survey confounding
# what is wrong with this?
fit0 <- sdmTMB(formula = dens ~ 0 + survey + fyear + gear + depth:survey,
               data=df, 
               time="year",
               spatial='off',
               spatiotemporal='off',
               mesh=mesh, 
               family = tweedie(link = "log"),
               control = sdmTMBcontrol(newton_loops = 1))
sanity(fit0)

fit0 <- sdmTMB(formula = dens ~ 0 + fyear + gear,
               data=df, 
               time="year",
               spatial='off',
               spatiotemporal='off',
               mesh=mesh, 
               family = tweedie(link = "log"),
               control = sdmTMBcontrol(newton_loops = 1))
sanity(fit0)

df$resids <- residuals(fit0)
png(file.path(plotdir, paste0('qqplot_fit0.png')), units='in', width=5, height=5, res=150)
qqnorm(df$resids)
abline(a = 0, b = 1)
dev.off()

tidy(fit0, "fixed", conf.int=T, digits=3) %>% print(n = Inf)
tidy(fit0, "ran_pars", conf.int=T, digits=3) %>% print(n = Inf)

# add spatial
fit1 <- sdmTMB(formula = dens ~ 0 + fyear + gear,
               data=df, 
               time="year",
               spatial='on',
               spatiotemporal='off',
               mesh=mesh, 
               family = tweedie(link = "log"),
               control = sdmTMBcontrol(newton_loops = 1))
sanity(fit1)

df$resids <- residuals(fit1)
png(file.path(plotdir, paste0('qqplot_fit1.png')), units='in', width=5, height=5, res=150)
qqnorm(df$resids)
abline(a = 0, b = 1)
dev.off()

tidy(fit1, "fixed", conf.int=T, digits=3) %>% print(n = Inf)
tidy(fit1, "ran_pars", conf.int=T, digits=3) %>% print(n = Inf)

# add spatiotemporal
fit2 <- sdmTMB(formula = dens ~ 0 + fyear + gear,
               data=df, 
               time="year",
               spatial='on',
               spatiotemporal='iid',
               mesh=mesh, 
               family = tweedie(link = "log"),
               control = sdmTMBcontrol(newton_loops = 1))
sanity(fit2)

df$resids <- residuals(fit2)
png(file.path(plotdir, paste0('qqplot_fit2.png')), units='in', width=5, height=5, res=150)
qqnorm(df$resids)
abline(a = 0, b = 1)
dev.off()

tidy(fit2, "fixed", conf.int=T, digits=3) %>% print(n = Inf)
tidy(fit2, "ran_pars", conf.int=T, digits=3) %>% print(n = Inf)

# delta gamma
fit3 <- sdmTMB(formula = dens ~ 0 + fyear + gear,
               data=df, 
               time="year",
               spatial='on',
               spatiotemporal='off',
               mesh=mesh, 
               family = delta_gamma(),
               control = sdmTMBcontrol(newton_loops = 1))
sanity(fit3)

df$resids <- residuals(fit3)
png(file.path(plotdir, paste0('qqplot_fit3.png')), units='in', width=5, height=5, res=150)
qqnorm(df$resids)
abline(a = 0, b = 1)
dev.off()

tidy(fit3, "fixed", conf.int=T, digits=3) %>% print(n = Inf)
tidy(fit3, "ran_pars", conf.int=T, digits=3) %>% print(n = Inf)

# delta gamma + depth
fit4 <- sdmTMB(formula = dens ~ 0 + fyear + gear + depth:survey,
               data=df, 
               time="year",
               spatial='on',
               spatiotemporal='off',
               mesh=mesh, 
               family = delta_gamma(),
               control = sdmTMBcontrol(newton_loops = 1))
sanity(fit4)

df$resids <- residuals(fit4)
png(file.path(plotdir, paste0('qqplot_fit4.png')), units='in', width=5, height=5, res=150)
qqnorm(df$resids)
abline(a = 0, b = 1)
dev.off()

tidy(fit4, "fixed", conf.int=T, digits=3) %>% print(n = Inf)
tidy(fit4, "ran_pars", conf.int=T, digits=3) %>% print(n = Inf)

# delta gamma + s(depth)
fit5 <- sdmTMB(formula = dens ~ 0 + fyear + gear + s(depth, by=survey),
               data=df, 
               time="year",
               spatial='on',
               spatiotemporal='off',
               mesh=mesh, 
               family = delta_gamma(),
               control = sdmTMBcontrol(newton_loops = 1))
sanity(fit5)

df$resids <- residuals(fit5)
png(file.path(plotdir, paste0('qqplot_fit5.png')), units='in', width=5, height=5, res=150)
qqnorm(df$resids)
abline(a = 0, b = 1)
dev.off()

tidy(fit5, "fixed", conf.int=T, digits=3) %>% print(n = Inf)
tidy(fit5, "ran_pars", conf.int=T, digits=3) %>% print(n = Inf)

# results from fit3
# spatial residuals
df$resids <- residuals(fit3)
png(file.path(plotdir, 'resids_spatial.png'), units='in', width=10, height=6, res=150)
print(ggplot(df, aes(X, Y, col = resids)) + scale_colour_gradient2() +
        geom_point() + facet_wrap(~year) + coord_fixed() + theme_bw())
dev.off()

# prediction grid
res <- 5 # in km
df$X_bin <- floor(df$X / res) * res
df$Y_bin <- floor(df$Y / res) * res
n_cells <- length(unique(paste0(df$X_bin, df$Y_bin)))
df$XY_bin <- paste(df$X_bin, df$Y_bin)
pred_grid <- data.frame(XY_bin = unique(df$XY_bin))
pred_grid$X <- as.numeric(unlist(lapply(strsplit(pred_grid$XY_bin," "), getElement, 1)))
pred_grid$Y <- as.numeric(unlist(lapply(strsplit(pred_grid$XY_bin," "), getElement, 2)))

# convert dataframe to sf object and set CRS 
pred_grid$X_m <- pred_grid$X*1000
pred_grid$Y_m <- pred_grid$Y*1000
df_sf <- st_as_sf(pred_grid, coords = c("X_m", "Y_m"), crs = get_crs(df, ll_names=c('lon','lat')))
df_geo <- st_transform(df_sf, crs = 4326) # transform to geographic coordinates (longitude/latitude)
coords <- st_coordinates(df_geo)

# Create a new dataframe with the coordinates
df_coords <- data.frame(
  lon = coords[, 1],
  lat = coords[, 2]
)
# check lat lon conversion
ggplot(df_coords, aes(lon, lat)) +
   geom_point() + theme_bw() 

# need all covariates, just use gear = Trawl
predgrid <- crossing(pred_grid, year=unique(df$year))
predgrid$gear <- 'Trawl'
predgrid$fyear <- factor(predgrid$year)
predgrid_geo <- crossing(df_coords, year=unique(df$year))
preds <- predict(fit3, newdata = predgrid, return_tmb_object = TRUE, type = "response")
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
index$year = as.integer(as.character(index$year))
g <- ggplot(index, aes(x=year, y=est)) + 
  ylab('Combined index (cod biomass)')   +
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

# compare to individual indices
index_trawl <- read.csv('exercises/coastal-survey-ex/plots/south/index_trawlsurvey_south.csv')
index_trawl$label = 'Trawl index'
index_garnruse <- readRDS('exercises/garn-ruse-ex/results/biomass-dg/inds_m2.rds')
index_garnruse$label = 'Garn ruse index'
index$label = 'Combined'

df <- rbind(index, index_trawl, index_garnruse)
df <- df %>% group_by(label) %>% mutate(est_stand = est / mean(est),
                                         lo_stand = lwr / mean(est),
                                         hi_stand = upr / mean(est))
years <- range(df$year)
png(file.path(plotdir,'compare_indices_standardized.png'), res=300, units='in', height=10, width=6)
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

png(file.path(plotdir,'compare_indices_standardized_nofacet.png'), res=300, units='in', height=6, width=7)
print(ggplot(df, aes(x=year, y=est_stand, ymin=lo_stand, ymax=hi_stand)) +
        geom_hline(yintercept=1, linetype=2) +
        geom_ribbon(aes(fill=label), alpha=0.15) +
        geom_line(aes(color=label), linewidth=1) + 
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
              legend.position = c(.85,.87)))
dev.off()

