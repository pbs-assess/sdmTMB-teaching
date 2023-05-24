# spatiotemporal index standardization of garn ruse survey
#   response: catch (number), ages independent, set level (not aggregated by site)

# sdmTMB vignette
# https://pbs-assess.github.io/sdmTMB/index.html

# 10-fold CV, block by site within year (site is the PSU, not set)
# run all models and compare log density
# Commander et al. 2022
# https://doi.org/10.7717/peerj.12783

# base: y ~ 0 + year + subarea + gear
m.formulas0 <- list(
  base = 'y ~ 0 + year + subarea', # no covariates or spatial field 
  base_spatial = 'y ~ 0 + year + subarea', # add spatial field
  base_st = 'y ~ 0 + year + subarea', # add spatiotemporal (iid) fields
  base_siteRE = 'y ~ 0 + year + subarea + (1 | site)', # RE of site instead of spatial field
  base_lat = 'y ~ 0 + year + subarea + lat', # linear effect of lat
  base_latsmooth = 'y ~ 0 + year + subarea + s(lat)', # spline on lat
  base_2dsmooth = 'y ~ 0 + year + subarea + s(lon,lat)' # 2d spline on lon,lat
)
m.names = names(m.formulas0)
m.formulas <- as.list(paste(m.formulas0, ' + gear'))
names(m.formulas) = m.names
n.mods = length(m.names)
m.s <- rep('off', n.mods) # default: no spatial field
m.s[grep('spatial', m.names)] = 'on' # add spatial field
m.s[grep('st', m.names)] = 'on' # add spatial field
m.st <- rep('off', n.mods) # default: no spatiotemporal fields
m.st[grep('st', m.names)] = 'iid' # add spatiotemporal (iid)

touse <- c("here","tidyverse","sdmTMB") 
lapply(touse, require, character.only=TRUE, quietly=TRUE)
source(here('exercises','garn-ruse-ex','helper-functions.R'))
plotdir <- here('exercises','garn-ruse-ex','plots')
resdir <- here('exercises','garn-ruse-ex','results')
if(!dir.exists(plotdir)) dir.create(plotdir, recursive=TRUE)
if(!dir.exists(resdir)) dir.create(resdir, recursive=TRUE)

# get clean data frame with catch-at-age
dfc <- readRDS(here('data','garn_ruse_cod_CAA.rds'))

# spatial model better on equal distance projection, not lat/lon
dfc <- add_utm_columns(dfc, c("lon", "lat"), units="km") # UTM zone 32 covers majority of data, 6-12E
mesh <- make_mesh(dfc, xy_cols = c("X","Y"), cutoff = 3)
png(file.path(plotdir, 'mesh_3km.png'), units='in', width=5, height=7, res=150)
plot(mesh)
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

# age 0 not caught in trammel/gillnets, don't fit effect of gear
df0 <- filter(dfc, gear == 'Fyke (ruse)')
mesh0 <- make_mesh(df0, xy_cols = c("X","Y"), cutoff = 3)

ages <- 0:5
n.ages <- length(ages)
years <- levels(dfc$year)
n.yr <- length(years)
fits <- vector('list', n.mods) # model fit objects
for(m in 1:n.mods){
  fits[[m]] <- vector('list', n.ages)
  for(a in 1:n.ages){
    # age 0 not caught in trammel/gillnets, don't fit effect of gear
    cat(paste0('model ',m,', age ',a))
    if(a == 1){
      df <- filter(dfc, gear == 'Fyke (ruse)')
      df$y = df[,paste0('catch_',a-1)] # total catch at age a
      fits[[m]][[a]] <- sdmTMB_cv(as.formula(m.formulas0[[m]]), 
                                  data=df, 
                                  time="year",
                                  fold_ids=df$fold,
                                  spatiotemporal=m.st[m], 
                                  spatial=m.s[m], 
                                  mesh=mesh0, 
                                  family = nbinom1())
    } else {
      df <- dfc
      df$y = df[,paste0('catch_',a-1)] # total catch at age a
      fits[[m]][[a]] <- sdmTMB_cv(as.formula(m.formulas[[m]]), 
                               data=df, 
                               time="year",
                               fold_ids=df$fold,
                               spatiotemporal=m.st[m], 
                               spatial=m.s[m], 
                               mesh=mesh, 
                               family = nbinom1())  
    }
    saveRDS(fits[[m]][[a]], file.path(resdir, paste0('fits_a',a,'_m',m,'.rds')))  
  }
}
# saveRDS(fits, file.path(resdir, 'fits.rds'))
# fits <- readRDS(file.path(resdir, 'fits.rds'))

res <- expand.grid(age = ages, model = m.names)
res[, c('pdHess.all','pdHess.n','maxgr.all','maxgr.n','elpd','loglik')] <- NA
for(m in 1:n.mods){
  for(a in 1:n.ages){
    fits <- readRDS(file.path(resdir, paste0('fits_a',a,'_m',m,'.rds')))
    ind = which(res$age == ages[a] & res$model == m.names[m])
    res$pdHess.all[ind] = fits$converged
    res$pdHess.n[ind] = sum(fits$pdHess)
    res$maxgr.all[ind] = all(fits$max_gradients < 0.001)
    res$maxgr.n[ind] = sum(fits$max_gradients < 0.001)
    res$elpd[ind] = fits$elpd
    res$loglik[ind] = fits$sum_loglik
  }
}

# compare models by elpd and loglik
# https://pbs-assess.github.io/sdmTMB/articles/web_only/cross-validation.html
# https://avehtari.github.io/modelselection/CV-FAQ.html#5_How_to_use_cross-validation_for_model_selection
# https://discourse.mc-stan.org/t/relationship-between-elpd-differences-and-likelihood-ratios/23855/4
# https://mc-stan.org/loo/reference/loo-glossary.html

elpd <- res %>% pivot_wider(id_cols=model, names_from=age, values_from=elpd) %>% as.data.frame
rownames(elpd) = elpd$model
elpd = elpd[,-1]
colnames(elpd) <- paste0('Age ',ages)
elpd$Total <- apply(elpd, 1, sum)
elpd.diff <- round(apply(elpd, 2, function(x) max(x) - x),3) # pos elpd, max is best
write.csv(elpd.diff, file=file.path(plotdir, 'elpd_diff.csv'))

elpd.diff
#                Age 0 Age 1 Age 2 Age 3 Age 4 Age 5 Total
# base           0.002 0.026 0.015 0.011 0.011 0.006 0.067
# base_spatial   0.005 0.002 0.004 0.003 0.001 0.002 0.014
# base_st        0.000 0.003 0.000 0.001 0.000 0.002 0.003
# base_siteRE    0.003 0.000 0.001 0.000 0.000 0.000 0.000
# base_lat       0.001 0.027 0.016 0.009 0.011 0.006 0.066
# base_latsmooth 0.001 0.016 0.014 0.006 0.012 0.004 0.050
# base_2dsmooth  0.004 0.011 0.011 0.002 0.011 0.004 0.040


loglik <- res %>% pivot_wider(id_cols=model, names_from=age, values_from=loglik) %>% as.data.frame
rownames(loglik) = loglik$model
loglik = loglik[,-1]
colnames(loglik) <- paste0('Age ',ages)
loglik$Total <- apply(loglik, 1, sum)
loglik.diff <- round(apply(loglik, 2, function(x) max(x) - x),3) # pos loglik, max is best
write.csv(loglik.diff, file=file.path(plotdir, 'loglik_diff.csv'))

loglik.diff
#                 Age 0  Age 1  Age 2  Age 3  Age 4  Age 5   Total
# base            0.749 79.326 12.445 17.763 24.491 11.280 130.439
# base_spatial    9.094  0.000  4.848  2.226  5.542  2.739   8.834
# base_st         3.509  4.328  5.880  5.451  6.896  2.739  13.187
# base_siteRE    10.280  1.909  0.000  0.000  0.000  3.426   0.000
# base_lat        0.000 81.049 13.810 15.786 25.626 11.060 131.715
# base_latsmooth  0.000 59.335 15.574 12.543 28.013  0.000  99.850
# base_2dsmooth   7.169 46.020 16.876  9.179 24.965  3.613  92.206

# ------------------------------------------------------------------------------
# fit full models
#   check diagnostics (qq plot, spatial residuals, )
#   compare difference in index
fits_full <- vector('list', n.mods) # model fit objects
preds <- vector('list', n.mods) # predicted catch at age for each year-site-gear combo, by model
inds <- vector('list', n.mods) # index by age, each model
res <- expand.grid(age = ages, model = m.names)
res[, c('pdHess','maxgr.ok','maxgr','no.na.se','aic')] <- NA
for(m in 1:n.mods){
  fits_full[[m]] <- vector('list', n.ages)
  for(a in 1:n.ages){
    # age 0 not caught in trammel/gillnets, don't fit effect of gear
    cat(paste0('model ',m,', age ',a))
    if(a == 1){
      df <- filter(dfc, gear == 'Fyke (ruse)')
      df$y = df[,paste0('catch_',a-1)] # total catch at age a
      fits_full[[m]][[a]] <- sdmTMB(as.formula(m.formulas0[[m]]), 
                                  data=df, 
                                  time="year",
                                  spatiotemporal=m.st[m], 
                                  spatial=m.s[m], 
                                  mesh=mesh0, 
                                  family = nbinom1())
    } else {
      df <- dfc
      df$y = df[,paste0('catch_',a-1)] # total catch at age a
      fits_full[[m]][[a]] <- sdmTMB(as.formula(m.formulas[[m]]), 
                                  data=df, 
                                  time="year",
                                  spatiotemporal=m.st[m], 
                                  spatial=m.s[m], 
                                  mesh=mesh, 
                                  family = nbinom1())  
    }
    saveRDS(fits_full[[m]][[a]], file.path(resdir, paste0('fits_full_a',a,'_m',m,'.rds')))  
    
    # convergence diagnostics
    tmp <- sanity(fits_full[[m]][[a]])
    save_sanity_text(tmp, fil=file.path(plotdir, paste0('convergence_age',a-1,'_m',m,'.txt')))
    
    tmp <- capture.output(fits_full[[m]][[a]]$sd_report)
    cat(tmp, sep="\n", file=file.path(plotdir, paste0('sdreport_age',a-1,'_m',m,'.txt')))
    
    # QQ plot residuals
    df$resids <- residuals(fits_full[[m]][[a]])
    png(file.path(plotdir, paste0('qqplot_age',a-1,'_m',m,'.png')), units='in', width=5, height=5, res=150)
    qqnorm(df$resids)
    abline(a = 0, b = 1)
    dev.off()
    
    # parameter estimates
    cat("Fixed effect parameters:
      ", sep="\n", file=file.path(plotdir, paste0('pars_age',a-1,'_m',m,'.txt')))
    tmp <- capture.output(tidy(fits_full[[m]][[a]], "fixed", conf.int=T))
    cat(tmp, sep="\n", file=file.path(plotdir, paste0('pars_age',a-1,'_m',m,'.txt')), append=TRUE)
    cat("
Random effect parameters:
      ", sep="\n", file=file.path(plotdir, paste0('pars_age',a-1,'_m',m,'.txt')), append=TRUE)
    tmp <- capture.output(tidy(fits_full[[m]][[a]], "ran_pars", conf.int=T))
    cat(tmp, sep="\n", file=file.path(plotdir, paste0('pars_age',a-1,'_m',m,'.txt')), append=TRUE)
    # tidy(fits_full[[a]], "ranef", conf.int=T)
    
    # marginal effects plots
    png(file.path(plotdir, paste0('effect_year_age',a-1,'_m',m,'.png')), units='in', width=5, height=5, res=150)
    print(visreg::visreg(fits_full[[m]][[a]], xvar = "year", scale='response', gg=T, rug=F, xlab='Year', ylab=paste0('Catch, age ',a-1)) + theme_bw())
    dev.off()
    
    png(file.path(plotdir, paste0('effect_subarea_age',a-1,'_m',m,'.png')), units='in', width=3, height=5, res=150)
    print(visreg::visreg(fits_full[[m]][[a]], xvar = "subarea", scale='response', gg=T, rug=F, xlab='Subarea', ylab=paste0('Catch, age ',a-1)) + theme_bw())
    dev.off()
    
    if(a > 1){
      png(file.path(plotdir, paste0('effect_gear_age',a-1,'_m',m,'.png')), units='in', width=5, height=5, res=150)
      print(visreg::visreg(fits_full[[m]][[a]], xvar = "gear", scale='response', gg=T, rug=F, xlab='Gear', ylab=paste0('Catch, age ',a-1)) + theme_bw())
      dev.off()
    }
    
    # spatial residuals
    png(file.path(plotdir, paste0('resids_spatial_age',a-1,'_m',m,'.png')), units='in', width=10, height=6, res=150)
    print(ggplot(df, aes(X, Y, col = resids)) + scale_colour_gradient2() +
            geom_point() + facet_wrap(~year) + coord_fixed() + theme_bw())
    dev.off()
    
    # prediction grid
    #   average lat/lon/depth per site and gear
    #   total (site) = 6*fyke + 2*trammel45 + 2*trammel37
    #   index = sum over sites
    predgrid <- df %>% group_by(site, gear) %>% 
      # summarize(X = mean(X), Y = mean(Y), subarea=unique(subarea))
      summarize(Hsig.ave=mean(Hsig.ave), Hbott.ave=mean(Hbott.ave), depth = mean(depth), X = mean(X), Y = mean(Y), subarea=unique(subarea), lat = mean(lat), lon = mean(lon))
    predgrid$area = 2
    predgrid$area[predgrid$gear == 'Fyke (ruse)'] = 6
    predgrid <- crossing(predgrid, year=unique(df$year))
    
    # predict at each year-site-gear combination (mean lat, lon, depth by site-gear)
    #   total (site) = 6*fyke + 2*trammel45 + 2*trammel37
    #   sum over sites to get index 
    preds[[m]][[a]] <- predict(fits_full[[m]][[a]], newdata = predgrid, return_tmb_object = TRUE)
    predgrid$est = preds[[m]][[a]]$data$est
    
    # plot predictions
    png(file.path(plotdir, paste0('preds_spatial_age',a-1,'_m',m,'.png')), units='in', width=10, height=6, res=150)
    print(plot_map(predgrid, "est") + scale_color_viridis_c() + ggtitle(paste0("Predicted catch age ",a-1)))
    dev.off()
    
    if(m.s[m] == 'on'){
      png(file.path(plotdir, paste0('spatial_field_age',a-1,'_m',m,'.png')), units='in', width=10, height=6, res=150)
      print(plot_map(preds[[m]][[a]]$data, "omega_s") + ggtitle(paste0("Spatial field age ",a-1)) + scale_color_gradient2() + guides(color = guide_colorbar(reverse = TRUE)))
      dev.off()
    }
    
    # index of abundance
    inds[[m]][[a]] <- get_index(preds[[m]][[a]], area=preds[[m]][[a]]$data$area, bias_correct = TRUE)
    inds[[m]][[a]]$year = as.integer(as.character(inds[[m]][[a]]$year))
    
    saveRDS(inds[[m]][[a]], file.path(resdir, paste0('inds_a',a-1,'_m',m,'.rds')))  
    saveRDS(preds[[m]][[a]], file.path(resdir, paste0('preds_a',a-1,'_m',m,'.rds')))
    
    # get results
    ind = which(res$age == ages[a] & res$model == m.names[m])
    res$pdHess[ind] = fits_full[[m]][[a]]$pos_def_hessian
    res$maxgr.ok[ind] = all(fits_full[[m]][[a]]$gradients < 0.001)
    res$maxgr[ind] =  max(fits_full[[m]][[a]]$gradients)
    res$no.na.se[ind] = all(!is.na(sqrt(diag(fits_full[[m]][[a]]$sd_report$cov.fixed))))
    res$aic[ind] = AIC(fits_full[[m]][[a]])
  }
}
saveRDS(res, file.path(resdir, 'res.rds'))
write.csv(res, file=file.path(plotdir, 'res.csv'), row.names = FALSE)
# res <- readRDS(file.path(resdir, 'res.rds'))

df.conv <- res %>% pivot_wider(id_cols=model, names_from=age, values_from=pdHess) %>% as.data.frame
rownames(df.conv) = df.conv$model
df.conv = df.conv[,-1]
# 0    1    2    3    4    5
# base           TRUE TRUE TRUE TRUE TRUE TRUE
# base_spatial   TRUE TRUE TRUE TRUE TRUE TRUE
# base_st        TRUE TRUE TRUE TRUE TRUE TRUE
# base_siteRE    TRUE TRUE TRUE TRUE TRUE TRUE
# base_lat       TRUE TRUE TRUE TRUE TRUE TRUE
# base_latsmooth TRUE TRUE TRUE TRUE TRUE TRUE
# base_2dsmooth  TRUE TRUE TRUE TRUE TRUE TRUE

df.aic <- res %>% pivot_wider(id_cols=model, names_from=age, values_from=aic) %>% as.data.frame
rownames(df.aic) = df.aic$model
df.aic = df.aic[,-1]
df.aic$total = apply(df.aic, 1, sum)
df.aic.diff <- round(apply(df.aic, 2, function(x) x - min(x)),3)
colnames(df.aic.diff) <- c(paste0('Age ',ages),'Total')
write.csv(df.aic.diff, file=file.path(plotdir, 'aic_diff.csv'))

#                 Age 0   Age 1  Age 2  Age 3  Age 4  Age 5   Total
# base           18.474 166.265 63.475 36.128 31.843 13.483 315.556
# base_spatial    5.989   2.765 13.863  3.048  7.979  0.791  20.322
# base_st         0.000   0.000  0.000  1.452  9.868  2.791   0.000
# base_siteRE     6.344   7.208 15.495  0.000  0.000  0.214  15.150
# base_lat       16.505 166.964 64.935 31.328 33.589 10.773 309.984
# base_latsmooth 18.505 134.511 63.522 26.248 35.589  0.151 264.415
# base_2dsmooth  18.018  80.109 57.716 20.264 31.280  0.000 193.276

# _st lowest AIC
# _siteRE highest elpd
# compare index and index CV of best 3
mods <- 2:4
ages <- 0:5
n.ages <- length(ages)
n.mods <- length(mods)
inds <- vector('list', n.mods)
ind.df <- data.frame(matrix(NA, nrow=n.yr*n.ages, ncol=6))
colnames(ind.df) <- c('year','est','lwr','upr','log_est','se')
for(m in 1:n.mods){
  for(a in 1:n.ages){
    inds[[m]][[a]] <- readRDS(file.path(resdir, paste0('inds_a',a-1,'_m',mods[m],'.rds')))
  }
  names(inds[[m]]) <- ages
}
names(inds) <- m.names[mods]
ind.df <- dplyr::bind_rows(lapply(inds, function(x) dplyr::bind_rows(x, .id = 'Age')), .id='Model')

png(file.path(plotdir, 'index_age_base_spatial3.png'), units='in', width=8, height=6, res=150)
print(plot_timeseries_compare(df=ind.df, addCI=TRUE, addCV=TRUE, 
                              val='est', facet.var='Age', color.var='Model',
                              relative=FALSE, ytitle='Index of abundance (numbers)', yzero=TRUE))
dev.off()

# For 'best' model plot
#  index at age, show larger CV in south years (fewer fish)
#  correlation adjacent year/age (consistency)
#  cohort tracking

# index at age, larger CV in south years (fewer fish)
ind.df2 <- filter(ind.df, Model == 'base_spatial')
ind.df2$Subarea = 'North'
ind.df2$Subarea[ind.df2$year %in% c(2015,2017,2019,2021)] = 'South'
ind.df2 %>% group_by(Subarea, Age) %>% 
  summarize(meanCV = mean(se)) %>% as.data.frame
png(file.path(plotdir, 'index_age_subareaCV.png'), units='in', width=8, height=6, res=150)
print(plot_timeseries_subareaCV(df=ind.df2))
dev.off()

# correlation adjacent year/age
plot_index_at_age_consistency_sdmtmb(ind.df2, col='est', plot=TRUE, fname=file.path(plotdir,'iaa-consistency.png'))
# same but on log scale
plot_index_at_age_consistency_sdmtmb(ind.df2, col='log_est', plot=TRUE, fname=file.path(plotdir,'iaa-consistency-log.png'))

# cohort tracking (johanna's plot)
agelabs <- paste('Age', unique(ind.df2$Age))
ind.df2$year = as.numeric(ind.df2$year)
col='est'
IAA <- ind.df2 %>% 
  select(Age, year, (!!as.symbol(col))) %>% 
  complete(year = min(year):max(year)) %>%
  pivot_wider(names_from = Age, values_from=(!!as.symbol(col))) %>%
  select(-c(year,`NA`)) %>%
  as.matrix
rownames(IAA) <- seq(min(years), max(years))
colnames(IAA) <- 0:5

png(file.path(plotdir,'cohorts-index0.png'), width = 18, height = 12, units = "cm",res=150)
plot.cohort.year(mat=IAA, age=0:5, start.coh=2011, stop.coh=max(ind.df2$year)-1, ylab='Index-at-age (numbers)')
dev.off()

png(file.path(plotdir,'cohorts-index1.png'), width = 18, height = 12, units = "cm",res=150)
plot.cohort.year(mat=IAA, age=1:5, start.coh=2011, stop.coh=max(ind.df2$year)-1, ylab='Index-at-age (numbers)')
dev.off()

# same but on log scale
col='log_est'
IAA <- ind.df2 %>% 
  select(Age, year, (!!as.symbol(col))) %>% 
  complete(year = min(year):max(year)) %>%
  pivot_wider(names_from = Age, values_from=(!!as.symbol(col))) %>%
  select(-c(year,`NA`)) %>%
  as.matrix
rownames(IAA) <- seq(min(years), max(years))
colnames(IAA) <- 0:5

png(file.path(plotdir,'cohorts-index0-log.png'), width = 18, height = 12, units = "cm",res=150)
plot.cohort.year(mat=IAA, age=0:5, start.coh=2011, stop.coh=max(ind.df2$year)-1, ylab='Log (Index-at-age, numbers)')
dev.off()

png(file.path(plotdir,'cohorts-index1-log.png'), width = 18, height = 12, units = "cm",res=150)
plot.cohort.year(mat=IAA, age=1:5, start.coh=2011, stop.coh=max(ind.df2$year)-1, ylab='Log (Index-at-age, numbers)')
dev.off()


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



