# spatiotemporal index standardization of garn ruse survey
# simpler example, get weight from L-W and then model total biomass per net, not number of fish at age

# still do manual block CV
# 10-fold CV, block by site within year (site is the PSU, not set)
# run all models and compare log density
# Commander et al. 2022
# https://doi.org/10.7717/peerj.12783

touse <- c("here","tidyverse","sdmTMB") 
lapply(touse, require, character.only=TRUE, quietly=TRUE)
source(here('exercises','garn-ruse-ex','helper-functions.R'))
plotdir <- here('exercises','garn-ruse-ex','plots','biomass')
resdir <- here('exercises','garn-ruse-ex','results','biomass')
if(!dir.exists(plotdir)) dir.create(plotdir, recursive=TRUE)
if(!dir.exists(resdir)) dir.create(resdir, recursive=TRUE)

# make base formula and then add different spatial models
base.formula <- 'biomass_kg ~ 0 + fyear + subarea + gear'
spatial.formulas <- list(
  base = '', # no covariates or spatial field 
  base_spatial = '', # add spatial field in sdmtmb function
  base_spatial_st = '', # add spatial + spatiotemporal (iid) fields in sdmtmb function
  base_siteRE = ' + (1 | site)', # RE of site instead of spatial field
  base_lat = ' + lat', # linear effect of lat
  base_latsmooth = ' + s(lat)', # spline on lat
  base_2dsmooth = ' + s(lon,lat)' # 2d spline on lon,lat
)
m.names = names(spatial.formulas)
m.formulas <- as.list(paste0(base.formula, spatial.formulas))
names(m.formulas) = m.names
m.formulas

n.mods = length(m.names)
m.s <- rep('off', n.mods) # default: no spatial field
m.s[grep('spatial', m.names)] = 'on' # add spatial field
m.st <- rep('off', n.mods) # default: no spatiotemporal fields
m.st[grep('st', m.names)] = 'iid' # add spatiotemporal (iid)

# get catch data frame, weight has been estimated from length using survey L-W
# each row = a net set, total cod biomass is 'biomass_kg'
dfc <- readRDS(here('data','garn_ruse_cod_biomass.rds'))
dfc$fyear <- factor(dfc$year)

# spatial model better on equal distance projection, not lat/lon
dfc <- add_utm_columns(dfc, c("lon", "lat"), units="km") # UTM zone 32 covers majority of data, 6-12E
mesh <- make_mesh(dfc, xy_cols = c("X","Y"), cutoff = 3)
png(file.path(plotdir, 'mesh_3km.png'), units='in', width=5, height=7, res=150)
plot(mesh)
dev.off()

# look at site structure of data
# table(dfc$site, dfc$year)
# zoom to 65-66 latitude
# notice:
#   1. sites (color) are consistent across years
#   2. not all sites sampled in all years
png(file.path(plotdir, 'site_year_65-66N.png'), units='in', width=10, height=7, res=150)
ggplot(filter(dfc, lat < 66 & lat > 65), aes(x=lon, y=lat, color=site)) + geom_point() + facet_wrap(vars(fyear)) + theme_bw()
dev.off()

# set up folds, 10 blocks by site within year
set.seed(1234)
foldfn <- function(x, nfolds){
  y = rep(NA, length(x))
  tmp <- split(sample(1:length(x)), rep(1:nfolds, length.out = length(x)))
  for(i in 1:nfolds) y[tmp[[i]]] = i
  return(y)
}
yrsite <- dfc %>% group_by(year, site) %>% distinct(site) %>% 
  ungroup() %>% group_by(year) %>% 
  mutate(fold = foldfn(site, nfolds=10)) %>% as.data.frame
dfc$fold = NA
for(i in 1:dim(dfc)[1]) dfc$fold[i] = yrsite$fold[yrsite$year == dfc$year[i] & yrsite$site == dfc$site[i]]

fits <- vector('list', n.mods) # model fit objects
for(m in 1:n.mods){
  cat(paste0('model ',m))
  fits[[m]] <- sdmTMB_cv(as.formula(m.formulas[[m]]), 
                               data=dfc, 
                               time="year",
                               fold_ids=dfc$fold,
                               spatiotemporal=m.st[m], 
                               spatial=m.s[m], 
                               mesh=mesh, 
                               family = tweedie(link = "log"))  
  saveRDS(fits[[m]], file.path(resdir, paste0('fits_m',m,'.rds')))
}
# saveRDS(fits, file.path(resdir, 'fits.rds'))
# fits <- readRDS(file.path(resdir, 'fits.rds'))

res <- data.frame(model = m.names)
res[, c('pdHess.all','pdHess.n','maxgr.all','maxgr.n','elpd','loglik')] <- NA
for(m in 1:n.mods){
    fit <- readRDS(file.path(resdir, paste0('fits_m',m,'.rds')))
    res$pdHess.all[m] = fit$converged
    res$pdHess.n[m] = sum(fit$pdHess) # number of folds
    res$maxgr.all[m] = all(fit$max_gradients < 0.001)
    res$maxgr.n[m] = sum(fit$max_gradients < 0.001) # number of folds
    res$elpd[m] = fit$elpd
    res$loglik[m] = fit$sum_loglik
}

# compare models by elpd and loglik
# https://pbs-assess.github.io/sdmTMB/articles/web_only/cross-validation.html
# https://avehtari.github.io/modelselection/CV-FAQ.html#5_How_to_use_cross-validation_for_model_selection
# https://discourse.mc-stan.org/t/relationship-between-elpd-differences-and-likelihood-ratios/23855/4
# https://mc-stan.org/loo/reference/loo-glossary.html

res$elpd.diff <- round(max(res$elpd) - res$elpd, 3) # max is best
res$loglik.diff <- round(max(res$loglik) - res$loglik, 3) # max is best
write.csv(res, file=file.path(plotdir, 'res_cv.csv'))

res
#             model pdHess.all pdHess.n maxgr.all maxgr.n      elpd    loglik elpd.diff loglik.diff
# 1            base       TRUE       10     FALSE       5 -1.527662 -3521.182     0.081      86.251
# 2    base_spatial       TRUE       10     FALSE       1 -1.448504 -3434.931     0.002       0.000
# 3 base_spatial_st       TRUE       10     FALSE       0 -1.446288 -3444.065     0.000       9.135
# 4     base_siteRE       TRUE       10     FALSE       1 -1.447648 -3439.472     0.001       4.541
# 5        base_lat       TRUE       10     FALSE       6 -1.521629 -3518.015     0.075      83.084
# 6  base_latsmooth       TRUE       10     FALSE       2 -1.506838 -3498.682     0.061      63.751
# 7   base_2dsmooth       TRUE       10     FALSE       7 -1.473998 -3469.855     0.028      34.924

# ------------------------------------------------------------------------------
# models 2-4 look close, fit to all data and compare
#   diagnostics (qq plot, spatial residuals, )
#   difference in index
mods = 2:4
n.mods <- length(mods)
fits_full <- vector('list', n.mods) # model fit objects
preds <- vector('list', n.mods) # predicted catch at age for each year-site-gear combo, by model
inds <- vector('list', n.mods) # index by age, each model
res <- data.frame(model = m.names[mods])
res[, c('pdHess','maxgr.ok','maxgr','no.na.se','aic')] <- NA
for(m in 1:n.mods){
  fits_full[[m]] <- sdmTMB(as.formula(m.formulas[[mods[m]]]), 
                         data=dfc, 
                         time="fyear",
                         spatiotemporal=m.st[mods[m]], 
                         spatial=m.s[mods[m]], 
                         mesh=mesh, 
                         family = tweedie(link = "log"),
                         control = sdmTMBcontrol(newton_loops = 1))  
  saveRDS(fits_full[[m]], file.path(resdir, paste0('fits_full_m',mods[m],'.rds')))
    
    # convergence diagnostics
    tmp <- sanity(fits_full[[m]])
    save_sanity_text(tmp, fil=file.path(plotdir, paste0('convergence_m',mods[m],'.txt')))
    
    tmp <- capture.output(fits_full[[m]]$sd_report)
    cat(tmp, sep="\n", file=file.path(plotdir, paste0('sdreport_m',mods[m],'.txt')))
    
    # QQ plot residuals
    dfc$resids <- residuals(fits_full[[m]])
    png(file.path(plotdir, paste0('qqplot_m',mods[m],'.png')), units='in', width=5, height=5, res=150)
    qqnorm(dfc$resids)
    abline(a = 0, b = 1)
    dev.off()
    
    # parameter estimates
    cat("Fixed effect parameters:
      ", sep="\n", file=file.path(plotdir, paste0('pars_m',mods[m],'.txt')))
    tmp <- capture.output(tidy(fits_full[[m]], "fixed", conf.int=T))
    cat(tmp, sep="\n", file=file.path(plotdir, paste0('pars_m',mods[m],'.txt')), append=TRUE)
    cat("
Random effect parameters:
      ", sep="\n", file=file.path(plotdir, paste0('pars_m',mods[m],'.txt')), append=TRUE)
    tmp <- capture.output(tidy(fits_full[[m]], "ran_pars", conf.int=T))
    cat(tmp, sep="\n", file=file.path(plotdir, paste0('pars_m',mods[m],'.txt')), append=TRUE)

    # marginal effects plots
    png(file.path(plotdir, paste0('effect_year_m',mods[m],'.png')), units='in', width=5, height=5, res=150)
    print(visreg::visreg(fits_full[[m]], xvar = "fyear", scale='response', gg=T, rug=F, xlab='Year', ylab='Biomass (kg)') + theme_bw())
    dev.off()
    
    png(file.path(plotdir, paste0('effect_subarea_m',mods[m],'.png')), units='in', width=3, height=5, res=150)
    print(visreg::visreg(fits_full[[m]], xvar = "subarea", scale='response', gg=T, rug=F, xlab='Subarea', ylab='Biomass (kg)') + theme_bw())
    dev.off()
    
    png(file.path(plotdir, paste0('effect_gear_m',mods[m],'.png')), units='in', width=5, height=5, res=150)
    print(visreg::visreg(fits_full[[m]], xvar = "gear", scale='response', gg=T, rug=F, xlab='Gear', ylab='Biomass (kg)') + theme_bw())
    dev.off()

    # spatial residuals
    png(file.path(plotdir, paste0('resids_spatial_m',mods[m],'.png')), units='in', width=10, height=6, res=150)
    print(ggplot(dfc, aes(X, Y, col = resids)) + scale_colour_gradient2() +
            geom_point() + facet_wrap(~fyear, nrow=2) + coord_fixed() + theme_bw())
    dev.off()
    
    # prediction grid
    #   average lat/lon/depth per site and gear
    #   total (site) = 6*fyke + 2*trammel45 + 2*trammel37
    #   index = sum over sites
    predgrid <- dfc %>% group_by(site, gear) %>% 
      # summarize(X = mean(X), Y = mean(Y), subarea=unique(subarea))
      summarize(depth = mean(depth), X = mean(X), Y = mean(Y), subarea=unique(subarea), lat = mean(lat), lon = mean(lon))
    predgrid$area = 2
    predgrid$area[predgrid$gear == 'Fyke (ruse)'] = 6
    predgrid <- crossing(predgrid, fyear=unique(dfc$fyear))
    
    # predict at each year-site-gear combination (mean lat, lon, depth by site-gear)
    #   total (site) = 6*fyke + 2*trammel45 + 2*trammel37
    #   sum over sites to get index 
    preds[[m]] <- predict(fits_full[[m]], newdata = predgrid, return_tmb_object = TRUE)
    predgrid$est = preds[[m]]$data$est
    
    # plot predictions
    png(file.path(plotdir, paste0('preds_spatial_m',mods[m],'.png')), units='in', width=10, height=6, res=150)
    print(plot_map_fyear(predgrid, "est") + scale_color_viridis_c() + ggtitle("Predicted biomass (kg)"))
    dev.off()
    
    if(m.s[mods[m]] == 'on'){
      png(file.path(plotdir, paste0('spatial_field_m',mods[m],'.png')), units='in', width=10, height=6, res=150)
      print(plot_map_fyear(preds[[m]]$data, "omega_s") + ggtitle("Spatial field") + scale_color_gradient2() + guides(color = guide_colorbar(reverse = TRUE)))
      dev.off()
    }
    
    # index of abundance
    inds[[m]] <- get_index(preds[[m]], area=preds[[m]]$data$area, bias_correct = TRUE)
    inds[[m]]$year = as.integer(as.character(inds[[m]]$fyear))
    
    saveRDS(inds[[m]], file.path(resdir, paste0('inds_m',mods[m],'.rds')))  
    saveRDS(preds[[m]], file.path(resdir, paste0('preds_m',mods[m],'.rds')))
    
    # get results
    res$pdHess[m] = fits_full[[m]]$pos_def_hessian
    res$maxgr.ok[m] = all(fits_full[[m]]$gradients < 0.001)
    res$maxgr[m] =  max(fits_full[[m]]$gradients)
    res$no.na.se[m] = all(!is.na(sqrt(diag(fits_full[[m]]$sd_report$cov.fixed))))
    res$aic[m] = AIC(fits_full[[m]])
}
res$daic <- round(res$aic - min(res$aic), 3)

saveRDS(res, file.path(resdir, 'res_fullfits.rds'))
write.csv(res, file=file.path(plotdir, 'res_fullfits.csv'), row.names = FALSE)

# _spatial highest CV loglik
# _st lowest AIC
# elpd diff very small
# compare index and index CV of best 3
names(inds) <- m.names[mods]
ind.df <- dplyr::bind_rows(inds, .id='Model')

png(file.path(plotdir, 'index_compare.png'), units='in', width=4, height=8, res=150)
print(plot_timeseries_compare(df=ind.df, addCI=TRUE, addCV=TRUE, 
                              val='est', facet.var='Model', color.var='Model',
                              relative=FALSE, ytitle='Biomass index', yzero=TRUE))
dev.off()
# 
# # For 'best' model plot
# #  index at age, show larger CV in south years (fewer fish)
# #  correlation adjacent year/age (consistency)
# #  cohort tracking
# 
# # index at age, larger CV in south years (fewer fish)
# ind.df2 <- filter(ind.df, Model == 'base_spatial')
# ind.df2$Subarea = 'North'
# ind.df2$Subarea[ind.df2$year %in% c(2015,2017,2019,2021)] = 'South'
# ind.df2 %>% group_by(Subarea, Age) %>% 
#   summarize(meanCV = mean(se)) %>% as.data.frame
# png(file.path(plotdir, 'index_age_subareaCV.png'), units='in', width=8, height=6, res=150)
# print(plot_timeseries_subareaCV(df=ind.df2))
# dev.off()
# 
# # correlation adjacent year/age
# plot_index_at_age_consistency_sdmtmb(ind.df2, col='est', plot=TRUE, fname=file.path(plotdir,'iaa-consistency.png'))
# # same but on log scale
# plot_index_at_age_consistency_sdmtmb(ind.df2, col='log_est', plot=TRUE, fname=file.path(plotdir,'iaa-consistency-log.png'))
# 
# # cohort tracking (johanna's plot)
# agelabs <- paste('Age', unique(ind.df2$Age))
# ind.df2$year = as.numeric(ind.df2$year)
# col='est'
# IAA <- ind.df2 %>% 
#   select(Age, year, (!!as.symbol(col))) %>% 
#   complete(year = min(year):max(year)) %>%
#   pivot_wider(names_from = Age, values_from=(!!as.symbol(col))) %>%
#   select(-c(year,`NA`)) %>%
#   as.matrix
# rownames(IAA) <- seq(min(years), max(years))
# colnames(IAA) <- 0:5
# 
# png(file.path(plotdir,'cohorts-index0.png'), width = 18, height = 12, units = "cm",res=150)
# plot.cohort.year(mat=IAA, age=0:5, start.coh=2011, stop.coh=max(ind.df2$year)-1, ylab='Index-at-age (numbers)')
# dev.off()
# 
# png(file.path(plotdir,'cohorts-index1.png'), width = 18, height = 12, units = "cm",res=150)
# plot.cohort.year(mat=IAA, age=1:5, start.coh=2011, stop.coh=max(ind.df2$year)-1, ylab='Index-at-age (numbers)')
# dev.off()
# 
# # same but on log scale
# col='log_est'
# IAA <- ind.df2 %>% 
#   select(Age, year, (!!as.symbol(col))) %>% 
#   complete(year = min(year):max(year)) %>%
#   pivot_wider(names_from = Age, values_from=(!!as.symbol(col))) %>%
#   select(-c(year,`NA`)) %>%
#   as.matrix
# rownames(IAA) <- seq(min(years), max(years))
# colnames(IAA) <- 0:5
# 
# png(file.path(plotdir,'cohorts-index0-log.png'), width = 18, height = 12, units = "cm",res=150)
# plot.cohort.year(mat=IAA, age=0:5, start.coh=2011, stop.coh=max(ind.df2$year)-1, ylab='Log (Index-at-age, numbers)')
# dev.off()
# 
# png(file.path(plotdir,'cohorts-index1-log.png'), width = 18, height = 12, units = "cm",res=150)
# plot.cohort.year(mat=IAA, age=1:5, start.coh=2011, stop.coh=max(ind.df2$year)-1, ylab='Log (Index-at-age, numbers)')
# dev.off()


# ------------------------------------------------------------------------------
# Extensions:

# # add covariates to 'spatial'
# depth = 'y ~ 0 + year + subarea + depth',
# Hbott = 'y ~ 0 + year + subarea + Hbott.ave',
# depth2 = 'y ~ 0 + year + subarea + depth + depth2',
# depth_Hbott = 'y ~ 0 + year + subarea + depth + Hbott.ave'

# note '+ gear' is added later because only for ages 1+
# m.formulas <- as.list(paste(m.formulas0, ' + gear'))

# estimate aggregate index as biomass of ages 2-5
# individual biological data, with weight interpolated using trawl survey LW
# df <- readRDS(here('data','garn_ruse_alk_expanded_LW.rds'))
# df$length_mm
# df$weight_g



