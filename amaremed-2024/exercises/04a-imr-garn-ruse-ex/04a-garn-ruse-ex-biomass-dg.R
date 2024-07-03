# spatiotemporal index standardization of garn ruse survey
# simpler example, get weight from L-W and then model total biomass per net, not number of fish at age

touse <- c("here","tidyverse","sdmTMB")
lapply(touse, require, character.only=TRUE, quietly=TRUE)
source(here('amaremed-2024','exercises','04a-imr-garn-ruse-ex','helper-functions.R'))
plotdir <- here('amaremed-2024','exercises','04a-imr-garn-ruse-ex','plots','biomass-dg')
resdir <- here('amaremed-2024','exercises','04a-imr-garn-ruse-ex','results','biomass-dg')
if(!dir.exists(plotdir)) dir.create(plotdir, recursive=TRUE)
if(!dir.exists(resdir)) dir.create(resdir, recursive=TRUE)

# make base formula and then add different spatial models
base.formula <- 'biomass_kg ~ 0 + fyear + subarea + gear + depth'
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
dfc <- readRDS(here('amaremed-2024','data','garn_ruse_cod_biomass.rds'))
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

fits_full <- vector('list', n.mods) # model fit objects
preds <- vector('list', n.mods) # predicted catch at age for each year-site-gear combo, by model
inds <- vector('list', n.mods) # index by age, each model
res <- data.frame(model = m.names)
res[, c('pdHess','maxgr.ok','maxgr','no.na.se','aic')] <- NA
for(m in 1:n.mods){
  fits_full[[m]] <- sdmTMB(as.formula(m.formulas[[m]]),
                         data=dfc,
                         time="year",
                         spatiotemporal=m.st[m],
                         spatial=m.s[m],
                         mesh=mesh,
                         family = delta_gamma(),
                         control = sdmTMBcontrol(newton_loops = 1))
  saveRDS(fits_full[[m]], file.path(resdir, paste0('fits_full_m',m,'.rds')))

    # convergence diagnostics
    tmp <- sanity(fits_full[[m]])
    save_sanity_text(tmp, fil=file.path(plotdir, paste0('convergence_m',m,'.txt')))

    tmp <- capture.output(fits_full[[m]]$sd_report)
    cat(tmp, sep="\n", file=file.path(plotdir, paste0('sdreport_m',m,'.txt')))

    # QQ plot residuals
    dfc$resids <- residuals(fits_full[[m]])
    png(file.path(plotdir, paste0('qqplot_m',m,'.png')), units='in', width=5, height=5, res=150)
    qqnorm(dfc$resids)
    abline(a = 0, b = 1)
    dev.off()

    # parameter estimates
    cat("Fixed effect parameters:
      ", sep="\n", file=file.path(plotdir, paste0('pars_m',m,'.txt')))
    tmp <- capture.output(tidy(fits_full[[m]], "fixed", conf.int=T))
    cat(tmp, sep="\n", file=file.path(plotdir, paste0('pars_m',m,'.txt')), append=TRUE)
    cat("
Random effect parameters:
      ", sep="\n", file=file.path(plotdir, paste0('pars_m',m,'.txt')), append=TRUE)
    tmp <- capture.output(tidy(fits_full[[m]], "ran_pars", conf.int=T))
    cat(tmp, sep="\n", file=file.path(plotdir, paste0('pars_m',m,'.txt')), append=TRUE)

    # marginal effects plots
    # png(file.path(plotdir, paste0('effect_year_m',m,'.png')), units='in', width=5, height=5, res=150)
    # print(visreg::visreg(fits_full[[m]], xvar = "fyear", scale='response', gg=T, rug=F, xlab='Year', ylab='Biomass (kg)') + theme_bw())
    # dev.off()
    #
    # png(file.path(plotdir, paste0('effect_subarea_m',m,'.png')), units='in', width=3, height=5, res=150)
    # print(visreg::visreg(fits_full[[m]], xvar = "subarea", scale='response', gg=T, rug=F, xlab='Subarea', ylab='Biomass (kg)') + theme_bw())
    # dev.off()
    #
    # png(file.path(plotdir, paste0('effect_gear_m',m,'.png')), units='in', width=5, height=5, res=150)
    # print(visreg::visreg(fits_full[[m]], xvar = "gear", scale='response', gg=T, rug=F, xlab='Gear', ylab='Biomass (kg)') + theme_bw())
    # dev.off()

    # spatial residuals
    png(file.path(plotdir, paste0('resids_spatial_m',m,'.png')), units='in', width=10, height=6, res=150)
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
    predgrid <- crossing(predgrid, year=unique(dfc$year))
    predgrid$fyear <- factor(predgrid$year)

    # predict at each year-site-gear combination (mean lat, lon, depth by site-gear)
    #   total (site) = 6*fyke + 2*trammel45 + 2*trammel37
    #   sum over sites to get index
    preds[[m]] <- predict(fits_full[[m]], newdata = predgrid, return_tmb_object = TRUE, type='response')
    predgrid$est = preds[[m]]$data$est

    # plot predictions
    png(file.path(plotdir, paste0('preds_spatial_m',m,'.png')), units='in', width=10, height=6, res=150)
    print(plot_map(predgrid, "est") + scale_color_viridis_c() + ggtitle("Predicted biomass (kg)"))
    dev.off()

    # if(m.s[m] == 'on'){
    #   png(file.path(plotdir, paste0('spatial_field_m',m,'.png')), units='in', width=10, height=6, res=150)
    #   print(plot_map(preds[[m]]$data, "omega_s") + ggtitle("Spatial field") + scale_color_gradient2() + guides(color = guide_colorbar(reverse = TRUE)))
    #   dev.off()
    # }

    # index of abundance
    inds[[m]] <- get_index(preds[[m]], area=preds[[m]]$data$area, bias_correct = TRUE)
    inds[[m]]$year = as.integer(as.character(inds[[m]]$year))

    saveRDS(inds[[m]], file.path(resdir, paste0('inds_m',m,'.rds')))
    saveRDS(preds[[m]], file.path(resdir, paste0('preds_m',m,'.rds')))

    # get results
    res$pdHess[m] = fits_full[[m]]$pos_def_hessian
    res$maxgr.ok[m] = all(fits_full[[m]]$gradients < 0.001)
    res$maxgr[m] =  max(fits_full[[m]]$gradients)
    res$no.na.se[m] = all(!is.na(sqrt(diag(fits_full[[m]]$sd_report$cov.fixed))))
    res$aic[m] = AIC(fits_full[[m]])
}
res$daic <- round(res$aic - min(res$aic), 3)
#             model pdHess maxgr.ok        maxgr no.na.se      aic    daic
# 1            base   TRUE     TRUE 1.317291e-07     TRUE 6736.676 187.104
# 2    base_spatial   TRUE     TRUE 2.052861e-07     TRUE 6551.715   2.143
# 3 base_spatial_st   TRUE     TRUE 5.709033e-08     TRUE 6552.441   2.869
# 4     base_siteRE   TRUE     TRUE 3.748516e-08     TRUE 6549.572   0.000
# 5        base_lat   TRUE     TRUE 4.445274e-04     TRUE 6714.323 164.751
# 6  base_latsmooth   TRUE     TRUE 2.396585e-07     TRUE 6691.760 142.188
# 7   base_2dsmooth   TRUE     TRUE 7.228533e-09     TRUE 6640.274  90.702

saveRDS(res, file.path(resdir, 'res_fullfits.rds'))
write.csv(res, file=file.path(plotdir, 'res_fullfits.csv'), row.names = FALSE)

# models 2-4 have similar AIC, compare index and index CV
mods <- 2:4
names(inds) <- m.names
ind.df <- dplyr::bind_rows(inds[mods], .id='Model')

png(file.path(plotdir, 'index_compare.png'), units='in', width=4, height=8, res=150)
print(plot_timeseries_compare(df=ind.df, addCI=TRUE, addCV=TRUE,
                              val='est', facet.var='Model', color.var='Model',
                              relative=FALSE, ytitle='Biomass index', yzero=TRUE))
dev.off()


