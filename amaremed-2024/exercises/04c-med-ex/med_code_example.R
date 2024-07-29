# AMAREMED 2024: Exercise 04c: Applications to Mediterranean data ----

# before install packages if needed
#install.packages(c("dplyr", "ggplot2","sdmTMB","maps","data.table","car","tidyr"))

library(dplyr)
library(ggplot2)
library(sdmTMB)
library(maps)
library(data.table)
library(car)
library(tidyr)

# Starting with hake example ----

# import and explore data
merl_TATBTC.AMMED <- readRDS("amaremed-2024/data/merl_TATBTC.AMMED.rds")
grid_1516 <- readRDS('amaremed-2024/data/grid_1516.rds')

mp <-map_data("world", region=c("Italy","Malta"))

#subset kg_km2: one value for each haul (code)
merl_kg <- merl_TATBTC.AMMED %>% distinct(kg_km2,code,area,codeIndex,X,Y,X.utm,Y.utm,year,month,PA,
                               depth,tmp_bot,sal_bot,tmp.sur,O2_bot,chl_int,poc_int,poc_bot,NO3_bot,PO4_bot)

# explore data
head(merl_kg)
dim(merl_kg)

g_kg_km2 <- ggplot() +
  geom_sf()+
  geom_point(merl_kg,mapping=aes(x=X,y=Y,size=kg_km2,color=as.factor(PA)))+
  geom_polygon(data=mp, aes(long,lat,group=group), fill="lightgray")+
  coord_sf(xlim = c(11, 15.5),
           ylim = c(35, 38.5),
           expand = TRUE)+theme_bw()+
  scale_x_continuous(breaks = c(11, 13,15))+
  scale_y_continuous(breaks = c(35,36, 37, 38))+ggtitle('Merluccius merluccius')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~year)
g_kg_km2

# Number of hauls in each year:
g.haul <- ggplot(merl_kg, aes(year))+ geom_bar()+theme_bw()+ggtitle('Number of hauls per year')
g.haul

# plot distribution of catches and look for outliers
g.kg_hist <- ggplot(merl_kg) +
  aes(x = kg_km2) + geom_histogram(
  bins = round(sqrt(length(merl_kg$kg_km2))), # set number of bins
  fill = "steelblue", color = "black")+ theme_minimal()
g.kg_hist

# histogram of log values
g.kg_hist_log <- ggplot(merl_kg) +
  aes(x = log(kg_km2+1)) + geom_histogram(
    bins = round(sqrt(length(merl_kg$kg_km2))), # set number of bins
    fill = "steelblue", color = "black")+ theme_minimal()
g.kg_hist_log

# TODO: decide whether we really need to remove outliers or not?
#1° try
out <- boxplot.stats(merl_kg$kg_km2)$out
out_ind <- which(merl_kg$kg_km2 %in% c(out))
merl_kg[out_ind,] # check outlier (above 0.75+1.5⋅IQR or below 025+1.5⋅IQR

#2° try
lower_bound <- quantile(merl_kg$kg_km2, 0.01)
upper_bound <- quantile(merl_kg$kg_km2, 0.99)
outlier_ind <- which(merl_kg$kg_km2 < lower_bound | merl_kg$kg_km2 > upper_bound)
merl_kg[outlier_ind, ]

#3° try
qqPlot(merl_kg$kg_km2) #  non normal distributed
qqPlot(log(merl_kg$kg_km2+1)) # non normal distributed
qqPlot(sqrt(merl_kg$kg_km2)) #  non normal distributed

# I try this solution
merl_kg <- merl_kg[-c(377,746),] #  better to change name

#env variable
matx.cor <- cor(merl_kg[,c('X','Y','depth','tmp_bot','sal_bot','tmp.sur','O2_bot','chl_int','poc_int','poc_bot','NO3_bot','PO4_bot','kg_km2')])
matx.cor2 <- ifelse(matx.cor<0.6 & matx.cor> -0.6, matx.cor ,NA)
#It is possible to change the variable in the model

# Normalize/scale predictor variables
merl_kg$depth <- scale(merl_kg$depth)
merl_kg$tmp_bot <- scale(merl_kg$tmp_bot)
merl_kg$sal_bot <- scale(merl_kg$sal_bot)
merl_kg$tmp.sur <- scale(merl_kg$tmp.sur)
merl_kg$O2_bot <- scale(merl_kg$O2_bot)
merl_kg$chl_int <- scale(merl_kg$chl_int)
merl_kg$poc_int <- scale(merl_kg$poc_int)
merl_kg$poc_bot <- scale(merl_kg$poc_bot)
merl_kg$NO3_bot <- scale(merl_kg$NO3_bot)
merl_kg$PO4_bot <- scale(merl_kg$PO4_bot)

#SDM
merl_kg$X.utm_km <- merl_kg$X.utm/1000
merl_kg$Y.utm_km <- merl_kg$Y.utm/1000
merl_kg$fyear <- as.factor(merl_kg$year)

mesh.merl_kg <- make_mesh(merl_kg, c("X.utm_km", "Y.utm_km"), cutoff = 10)
plot(mesh.merl_kg)

# compare models going from high to low complexity of spatiotemporal structure ----
# (as an idditional exercise, you can do the same comparing models with different variables)

# first, most complex model takes several minutes to fit
# (consider reducing mesh complexity to increase speed)
fit.data_kgkm2 <- sdmTMB(
  kg_km2 ~ 0+s(tmp_bot)+s(poc_bot)+s(depth)+fyear,
  #spatial_varying = ~ 0 + scaled_year, # I tried to apply the temporal trend, but it's very expensive from the calculation time
  family = tweedie(link = "log"),
  data = merl_kg,
  mesh = mesh.merl_kg,
  time = "year",
  spatial = "on",
  spatiotemporal = "ar1"
)
#saveRDS(fit.data_kgkm2,'amaremed-2024/exercises/04c-med-ex/fit.data_kgkm2.rds')

fit2.data_kgkm2 <- sdmTMB(
  kg_km2 ~ 0+s(tmp_bot)+s(poc_bot)+s(depth)+fyear,
  family = tweedie(link = "log"),
  data = merl_kg,
  mesh = mesh.merl_kg,
  time = "year",
  spatial = "on",
  spatiotemporal = "iid"
)
#saveRDS(fit.data_kgkm2,'amaremed-2024/exercises/04c-med-ex/fit2.data_kgkm2.rds')

fit3.data_kgkm2 <- sdmTMB(
  kg_km2 ~ 0+s(tmp_bot)+s(poc_bot)+s(depth)+fyear,
  family = tweedie(link = "log"),
  data = merl_kg,
  mesh = mesh.merl_kg,
  time = "year",
  spatial = "on",
  spatiotemporal = "off"
)
#saveRDS(fit.data_kgkm2,'amaremed-2024/exercises/04c-med-ex/fit3.data_kgkm2.rds')

fit4.data_kgkm2 <- sdmTMB(
  kg_km2 ~ 0+s(tmp_bot)+s(poc_bot)+s(depth)+fyear,
  family = tweedie(link = "log"),
  data = merl_kg,
  mesh = mesh.merl_kg,
  time = "year",
  spatial = "off",
  spatiotemporal = "off"
)
#saveRDS(fit.data_kgkm2,'amaremed-2024/exercises/04c-med-ex/fit4.data_kgkm2.rds')

# Compare fit and model performance, AIC across models ----
# Marginal AIC is not as good as cross-validation for model comparison, but here we use it for speed
AIC(fit.data_kgkm2)
AIC(fit2.data_kgkm2)
AIC(fit3.data_kgkm2)
AIC(fit4.data_kgkm2)

# examples of other fit metrics applied to model 1 above
# Exercise: try replicating these and other residual diagnostics for the other models too
pred.data <- predict(fit.data_kgkm2,type='response',return_tmb_object = TRUE)
plot(pred.data$data$kg_km2, pred.data$data$est)
rsq <- function (x, y) cor(x, y) ^ 2
rsq(pred.data$data$kg_km2, pred.data$data$est)


# PREDICTING AND COMPUTING DERIVED INDICES, using example of model 1 above ----
# Visualizing prediction on a continuous grid ----
grid_1516$X.utm_km <- grid_1516$X.utm/1000
grid_1516$Y.utm_km <- grid_1516$Y.utm/1000
grid_1516$fyear <- as.factor(grid_1516$year)

# Normalize/scale predictor variables on grid
grid_1516$depth <- scale(grid_1516$depth)
grid_1516$tmp_bot <- scale(grid_1516$tmp_bot)
grid_1516$sal_bot <- scale(grid_1516$sal_bot)
grid_1516$tmp.sur <- scale(grid_1516$tmp.sur)
grid_1516$O2_bot <- scale(grid_1516$O2_bot)
grid_1516$chl_int <- scale(grid_1516$chl_int)
grid_1516$poc_int <- scale(grid_1516$poc_int)
grid_1516$poc_bot <- scale(grid_1516$poc_bot)
grid_1516$NO3_bot <- scale(grid_1516$NO3_bot)
grid_1516$PO4_bot <- scale(grid_1516$PO4_bot)

pred.grid <- predict(fit.data_kgkm2,grid_1516[grid_1516$year<=2021,],type='response',return_tmb_object = TRUE)
max.est_kg <- round(max(pred.grid$data$est),2)
#saveRDS(pred.grid,'amaremed-2024/exercises/04c-med-ex/pred.grid_kg.rds')

g.map_merl_kg <- ggplot() +
  geom_sf()+
  geom_polygon(data=mp, mapping=aes(long,lat,group=group), fill="lightgray")+
  #geom_tile(pred.grid$data,mapping=aes(x=X,y=Y,fill=est))+
  geom_raster(pred.grid$data, mapping=aes(x=X,y=Y,fill=est))+
  scale_fill_viridis_c()+
  coord_sf(xlim = c(11, 15.5),
           ylim = c(35, 38.5),
           expand = TRUE)+theme_bw()+
  scale_x_continuous(breaks = c(11, 13,15))+
  scale_y_continuous(breaks = c(35,36, 37, 38))+
  ggtitle(paste('Merluccius merluccius','max estimation=', max.est_kg,'kg/km^2', sep=' '))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~year)
g.map_merl_kg
ggsave("amaremed-2024/exercises/04c-med-ex/plots/g.map_merl_kg.png", plot=g.map_merl_kg, width=7,height=7)

# get and plot index of abundance
ind.kg_merl <- get_index(pred.grid, level = 0.8, # to get 80% CI
                         area = 1, # could include actual area of grid cells for abs abundance estimate
                         bias_correct = FALSE) # bias correction off for speed
p.ind.kg <- ggplot(ind.kg_merl, aes(year, est)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab('Year') + ylab('Biomass estimate (kg)')
p.ind.kg
ggsave("amaremed-2024/exercises/04c-med-ex/plots/p.ind.kg.png", plot=p.ind.kg, width=7,height=7)

# get and plot center of gravity in km eastings and northings
cog.kg_merl <- get_cog(pred.grid, level = 0.8, area = 1, format = "wide")
p.cog.kg <- ggplot(cog.kg_merl, aes(est_x, est_y, colour = year)) +
  geom_point() +
  geom_linerange(aes(xmin = lwr_x, xmax = upr_x)) +
  geom_linerange(aes(ymin = lwr_y, ymax = upr_y)) +
  scale_colour_viridis_c()+ggtitle('Center of Gravity 1999:2021')+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   text=element_text(size=16))
p.cog.kg
ggsave("amaremed-2024/exercises/04c-med-ex/plots/p.cog.kg.png", plot=p.cog.kg, width=7,height=7)

#EXERCISE: fit models with month as a covariate to account for seasonal change in sampling
# options include modeling this as an factor
# or perhaps a cyclic smoother bs = 'cc'


# Shrimp example ----

# import and explore data
pape_TATBTC.AMMED <- readRDS('amaremed-2024/data/pape_TATBTC.AMMED.rds')

#subset kg_km2: one value for each haul (code)
pape_kg <- merl_TATBTC.AMMED %>% distinct(kg_km2,code,area,codeIndex,X,Y,X.utm,Y.utm,year,month,PA,
                                          depth,tmp_bot,sal_bot,tmp.sur,O2_bot,chl_int,poc_int,poc_bot,NO3_bot,PO4_bot)

head(pape_kg)
dim(pape_kg)

g_kg_km2 <- ggplot() +
  geom_sf()+
  geom_point(pape_kg,mapping=aes(x=X,y=Y,size=kg_km2,color=as.factor(PA)))+
  geom_polygon(data=mp, aes(long,lat,group=group), fill="lightgray")+
  coord_sf(xlim = c(11, 15.5),
           ylim = c(35, 38.5),
           expand = TRUE)+theme_bw()+
  scale_x_continuous(breaks = c(11, 13,15))+
  scale_y_continuous(breaks = c(35,36, 37, 38))+ggtitle('Parapenaeus longirostris')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~year)
g_kg_km2

# Number of hauls in each year:
g.haul <- ggplot(pape_kg, aes(year))+ geom_bar()+theme_bw()+ggtitle('Number of hauls per year')
g.haul

# plot distribution of catches and look for outliers
g.kg_hist <- ggplot(pape_kg) +
  aes(x = kg_km2) + geom_histogram(
    bins = round(sqrt(length(pape_kg$kg_km2))), # set number of bins
    fill = "steelblue", color = "black")+ theme_minimal()
g.kg_hist

# histogram of log values
g.kg_hist_log <- ggplot(pape_kg) +
  aes(x = log(kg_km2+1)) + geom_histogram(
    bins = round(sqrt(length(pape_kg$kg_km2))), # set number of bins
    fill = "steelblue", color = "black")+ theme_minimal()
g.kg_hist_log

# data prep
pape_kg$X.utm_km <- pape_kg$X.utm/1000
pape_kg$Y.utm_km <- pape_kg$Y.utm/1000
pape_kg$fyear <- as.factor(pape_kg$year)

mesh.pape_kg <- make_mesh(pape_kg, c("X.utm_km", "Y.utm_km"), cutoff = 10)
plot(mesh.pape_kg)

# fit model for shrimp
# EXERCISE: fit alternative models with different covariates (after scaling) and
# spatiotemporal settings and compare them
fit5.data_kgkm2 <- sdmTMB(
  kg_km2 ~ 0+fyear,
  family = tweedie(link = "log"),
  data = pape_kg,
  mesh = mesh.pape_kg,
  time = "year",
  spatial = "on",
  spatiotemporal = "ar1"
)
#saveRDS(fit5.data_kgkm2,'amaremed-2024/exercises/04c-med-ex/fit5.data_kgkm2.rds')

# predict to grid
pred.grid.pape <- predict(fit5.data_kgkm2,grid_1516[grid_1516$year<=2021,],type='response',return_tmb_object = TRUE)
#saveRDS(pred.grid.pape,'amaremed-2024/exercises/04c-med-ex/pred.grid.pape_kg.rds')

g.map_pape_kg <- ggplot() +
  geom_sf()+
  geom_polygon(data=mp, mapping=aes(long,lat,group=group), fill="lightgray")+
  #geom_tile(pred.grid$data,mapping=aes(x=X,y=Y,fill=est))+
  geom_raster(pred.grid$data, mapping=aes(x=X,y=Y,fill=est))+
  scale_fill_viridis_c()+
  coord_sf(xlim = c(11, 15.5),
           ylim = c(35, 38.5),
           expand = TRUE)+theme_bw()+
  scale_x_continuous(breaks = c(11, 13,15))+
  scale_y_continuous(breaks = c(35,36, 37, 38))+
  ggtitle(paste('Parapenaeus longirostris','max estimation=', max.est_kg,'kg/km^2', sep=' '))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~year)
g.map_pape_kg
ggsave("amaremed-2024/exercises/04c-med-ex/plots/g.map_pape_kg.png", plot=g.map_pape_kg, width=7,height=7)


# Calculate spatial overlap between hake and shrimp over time, here using the
# Bhattacharyya coefficient of the "closeness" of two random statistical samples
p1 <- pred.grid$data %>%
  rename(est_merl = est)
p2 <- pred.grid.pape$data %>%
  rename(est_pape = est)

pjoint <- p1 %>%
  left_join(p2,by=c("year","index")) %>%
  group_by(year) %>%
  mutate(bhat=sum(sqrt(est_merl*est_pape))) %>%
  summarize(mean_bhat = mean(bhat, na.rm = TRUE))
pjoint

## NOTE: I may have butchered this particular overlap metric but you get the gist,
## also, I suspect these values may be very high because of the difference in
## the scale of abundance between the hake and shrimp. There are probably easy
## ways to standardize that for the comparison (e.g., at the coarsest scale you
## might compare the predictions for each species based on overlap in their probability
## of occurrence or some other habitat suitability metric
#jpeg("amaremed-2024/exercises/04c-med-ex/plots/overlap_index.png")
plot(pjoint$year, pjoint$mean_bhat, type = "l",
     pch = 16, frame = FALSE,
     xlab = "Year", ylab = "Spatial overlap of biomass", col = "blue")
#dev.off()


# EXERCISE: calculate and plot index of abundance and COG for shrimp

# BONUS MULTISPECIES EXERCISE:
# try fitting models where the other species is a covariate
# e.g., hake_kg ~ shrimp_kg + ...





# Extra code to split models by ontogenetic stage ----
# 'Merl N/km2'
# #based on literature references: MEDISEH & GFCM_SAF_2021
#
# 'R	<15	cm	MERLMER' # recruits
# 'S	>=33	cm	MERLMER' # spawners
# 'J	>=15 <33	cm	MERLMER' # juveniles
#
# merl_TATBTC.AMMED$stage <- ifelse(merl_TATBTC.AMMED$length_class<150,'R', # length_class is in mm in the data set
#                             ifelse(merl_TATBTC.AMMED$length_class>=150 & merl_TATBTC.AMMED$length_class <330,'J',
#                              ifelse(merl_TATBTC.AMMED$length_class>=330,'S',NA)))
#
# # execute this code again, because before we eliminate the outlier
# merl_kg <- merl_TATBTC.AMMED %>% distinct(kg_km2,code,area,codeIndex,X,Y,X.utm,Y.utm,year,month,PA,
#             depth,tmp_bot,sal_bot,tmp.sur,O2_bot,chl_int,poc_int,poc_bot,NO3_bot,PO4_bot)
#
#
# unique(merl_TATBTC.AMMED$stage)
# agg.N_km2 <- aggregate( N_km2 ~ X+Y+code+stage,FUN=sum,data=merl_TATBTC.AMMED) # aggregate the number of individuals for stage (J,R or S)
#
# # each haul contains 0 or one of the different levels or a combination of the three levels J,R,S
# expand.merl <- agg.N_km2 %>%
#   expand(code,stage) %>%
#   merge(merl_kg,all=T,by='code') # use merl_kg to merge all variable included
#
# merge_merl <- agg.N_km2 %>% dplyr::right_join(expand.merl)
#
# merge_merl$X.utm <- merge_merl$X.utm/1000
# merge_merl$Y.utm <- merge_merl$Y.utm/1000
#
# #split data set for stage
# spl.merl <- split(merge_merl,merge_merl$stage)
# dim(spl.merl$J)
# dim(spl.merl$R)
# dim(spl.merl$S)
#
# merl_juv <- spl.merl$J #juveniles
# merl_rec <- spl.merl$R # recruits
# merl_spaw <- spl.merl$S # spawners
#
# #Before the following part, we can build a superfunction that applies to any size,
# #but I've kept the code reliable for every single data set, it's up to you!
#
#
# 'juvenile'
# #3° try
# qqPlot(merl_juv$N_km2) #  non normal distributed
# qqPlot(log(merl_juv$N_km2+1)) # non normal distributed
# qqPlot(sqrt(merl_juv$N_km2)) #  non normal distributed
#
# # I try this solution
# merl_juv <- merl_juv[-c(289,1574),] #  better to change name
#
# mesh.juv <- make_mesh(merl_juv, c("X.utm", "Y.utm"),cutoff = 2)
# plot(mesh.juv)
#
# fit_merl.juv <- sdmTMB(
#   N_km2 ~ s(tmp_bot, k = 5)+s(sal_bot, k = 5)+s(O2_bot,k=5)+s(poc_bot,k=5)+s(NO3_bot,k=5),
#   data = merl_juv,
#   mesh = mesh.juv,
#   time = "year",
#   family = tweedie(link = "log"),
#   spatial = "off",
#   spatiotemporal = "ar1"
# )
#
# saveRDS(fit_merl.juv,'fit_merl.juv.rds')
#
# # qqplot residuals
# rq_juv <- residuals(fit_merl.juv,type = 'mle-eb') #mle-mvn
# rq_juv <- rq_juv[is.finite(rq_juv)] # some Inf
# qqnorm(rq_juv);qqline(rq_juv)
#
# 'diagnostic juv'
# pred.data.juv <- predict(fit_merl.juv,type='response',return_tmb_object = TRUE)
# plot(pred.data.juv$data$N_km2, pred.data.juv$data$est)
# rsq <- function (x, y) cor(x, y) ^ 2
# rsq(pred.data.juv$data$N_km2, pred.data.juv$data$est)
# AIC(fit_merl.juv) #
#
# 'pred grid juv'
# pred.grid.juv <- predict(fit_merl.juv,grid_1516[grid_1516$year<=2021,],type='response',return_tmb_object = TRUE)
# max.est_juv <- round(max(pred.grid.juv$data$est),2)
# saveRDS(pred.grid.juv,'pred.grid.juv.rds')
#
# g.map_merl_juv <- ggplot() +
#   geom_sf()+
#   geom_polygon(data=mp, mapping=aes(long,lat,group=group), fill="lightgray")+
#   #geom_tile(pred.grid$data,mapping=aes(x=X,y=Y,fill=est))+
#   geom_raster(pred.grid.juv$data, mapping=aes(x=X,y=Y,fill=est))+
#   #geom_point(data=merl_kg,mapping=aes(x=X,y=Y,size=kg_km2),shape=21,color='red')+ too much!
#   scale_fill_viridis_c(limit=c(0,2000))+
#   coord_sf(xlim = c(11, 15.5),
#            ylim = c(35, 38.5),
#            expand = TRUE)+theme_bw()+
#   scale_x_continuous(breaks = c(11, 13,15))+
#   scale_y_continuous(breaks = c(35,36, 37, 38))+ggtitle(paste('Merluccius merluccius','max estimation=', max.est_juv,'N/km^2', sep=' '))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   facet_wrap(~year)
# ggsave("g.map_merl_juv.png", plot=g.map_merl_juv, width=7,height=7)
#
#
#
# 'recruits'
# #3° try
# qqPlot(merl_rec$N_km2) #  non normal distributed
# qqPlot(log(merl_rec$N_km2+1)) # non normal distributed
# qqPlot(sqrt(merl_rec$N_km2)) #  non normal distributed
#
# # I try this solution
# merl_rec <- merl_rec[-c(398,143),] #  better to change name
#
# mesh.rec <- make_mesh(merl_rec, c("X.utm", "Y.utm"),cutoff = 2)
# plot(mesh.rec)
#
# fit_merl.rec <- sdmTMB(
#   N_km2 ~ s(tmp_bot, k = 5)+s(sal_bot, k = 5)+s(O2_bot,k=5)+s(poc_bot,k=5)+s(NO3_bot,k=5),
#   data = merl_rec,
#   mesh = mesh.rec,
#   time = "year",
#   family = tweedie(link = "log"),
#   spatial = "off",
#   spatiotemporal = "ar1"
# )
#
# 'In this case we have an error in the fitting of the recruits'
#
# 'diagnostic rec'
# pred.data.rec <- predict(fit_merl.rec,type='response',return_tmb_object = TRUE)
# plot(pred.data.rec$data$N_km2, pred.data.rec$data$est)
# rsq <- function (x, y) cor(x, y) ^ 2
# rsq(pred.data.rec$data$N_km2, pred.data.rec$data$est)
# AIC(fit_merl.rec) #
#
# 'pred grid rec'
# pred.grid.rec <- predict(fit_merl.rec,grid_1516[grid_1516$year<=2021,],type='response',return_tmb_object = TRUE)
# max.est_rec <- round(max(pred.grid.rec$data$est),2)
# saveRDS(pred.grid.rec,'pred.grid.rec.rds')
#
# g.map_merl_rec <- ggplot() +
#   geom_sf()+
#   geom_polygon(data=mp, mapping=aes(long,lat,group=group), fill="lightgray")+
#   #geom_tile(pred.grid$data,mapping=aes(x=X,y=Y,fill=est))+
#   geom_raster(pred.grid.rec$data, mapping=aes(x=X,y=Y,fill=est))+
#   #geom_point(data=merl_kg,mapping=aes(x=X,y=Y,size=kg_km2),shape=21,color='red')+ too much!
#   scale_fill_viridis_c(limit=c(0,2000))+
#   coord_sf(xlim = c(11, 15.5),
#            ylim = c(35, 38.5),
#            expand = TRUE)+theme_bw()+
#   scale_x_continuous(breaks = c(11, 13,15))+
#   scale_y_continuous(breaks = c(35,36, 37, 38))+ggtitle(paste('Merluccius merluccius','max estimation=', max.est_rec,'N/km^2', sep=' '))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   facet_wrap(~year)
# ggsave("g.map_merl_rec.png", plot=g.map_merl_rec, width=7,height=7)
#
#
#
# 'and so on for spawners......'
#
# 'spawners'
# #3° try
# qqPlot(merl_spaw$N_km2) #  non normal distributed
# qqPlot(log(merl_spaw$N_km2+1)) # non normal distributed
# qqPlot(sqrt(merl_spaw$N_km2)) #  non normal distributed
#
# # I try this solution
# merl_spaw <- merl_spaw[-c(184,158),] #  better to change name
#
# mesh.spaw <- make_mesh(merl_spaw, c("X.utm", "Y.utm"),cutoff = 2)
# plot(mesh.spaw)
#
# fit_merl.spaw <- sdmTMB(
#   N_km2 ~ s(tmp_bot, k = 5)+s(sal_bot, k = 5)+s(O2_bot,k=5)+s(poc_bot,k=5)+s(NO3_bot,k=5),
#   data = merl_spaw,
#   mesh = mesh.spaw,
#   time = "year",
#   family = tweedie(link = "log"),
#   spatial = "off",
#   spatiotemporal = "ar1"
# )
#
# saveRDS(fit_merl.spaw,'fit_merl.spaw.rds')
#
#
#
# 'pred grid spaw'
# pred.grid.spaw <- predict(fit_merl.spaw,grid_1516[grid_1516$year<=2021,],type='response',return_tmb_object = TRUE)
# max.est_spaw <- round(max(pred.grid.spaw$data$est),2)
# saveRDS(pred.grid.spaw,'pred.grid.spaw.rds')
#
# g.map_merl_spaw <- ggplot() +
#   geom_sf()+
#   geom_polygon(data=mp, mapping=aes(long,lat,group=group), fill="lightgray")+
#   #geom_tile(pred.grid$data,mapping=aes(x=X,y=Y,fill=est))+
#   geom_raster(pred.grid.spaw$data, mapping=aes(x=X,y=Y,fill=est))+
#   #geom_point(data=merl_kg,mapping=aes(x=X,y=Y,size=kg_km2),shape=21,color='red')+ too much!
#   scale_fill_viridis_c(limit=c(0,80))+
#   coord_sf(xlim = c(11, 15.5),
#            ylim = c(35, 38.5),
#            expand = TRUE)+theme_bw()+
#   scale_x_continuous(breaks = c(11, 13,15))+
#   scale_y_continuous(breaks = c(35,36, 37, 38))+ggtitle(paste('Merluccius merluccius','max estimation=', max.est_spaw,'N/km^2', sep=' '))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   facet_wrap(~year)
# ggsave("g.map_merl_spaw.png", plot=g.map_merl_spaw, width=7,height=7)
#
