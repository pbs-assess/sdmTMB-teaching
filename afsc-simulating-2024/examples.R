library(sdmTMB)
library(ggplot2)

set.seed(123)

# Simulating from scratch ----

# scale number of primary sampling units in sampling universe
n <- 1000 # produces 10,000 per year
# set number of years
nyears <- 6
  
# make fake predictor (depth) and sampling locations:
predictor_dat <- data.frame(
  expand.grid(X = 1:(n/10), Y = 1:(n/10)),
  depth = rep(c(1:(n/20), rev(1:(n/20))), n/10) * rep(c(1,1.1,0.9,1.2,0.95,0.8), each = n*10), # * rlnorm(n*10, sdlog = 0.1) 
  year = rep(1:nyears, each = n*10)            
)
#saveRDS(predictor_dat, "grid_depth.RDS")

# plot domain for single year
ggplot(dplyr::filter(predictor_dat, year == 1), aes(X, Y)) +
  geom_tile(aes(fill = depth)) +
  scale_color_gradient2()
#ggsave("depth.png")

#---- fit model
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), type = "cutoff_search", n_knots = 100)

ranges <- seq(25, 100, 25) # specify values of range to simulate from

for(i in 1:length(ranges)){
  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + depth,
    data = predictor_dat,
    mesh = mesh,
    family = gaussian(link = "identity"),
    time = "year",
    B = c(0.1, 0.7), # B0 = intercept, B1 = depth coefficient slope
    range = ranges[i],
    rho = NULL, # or value between -1 and 1 to add spatiotemporal correlation
    sigma_O = 0.2,
    sigma_E = 10,
    phi = 5, # SD of observation error in this (Gaussian) case
    seed = 42
  )

  sim_dat$eta <- (sim_dat$eta - min(sim_dat$eta))
  sim_dat$observed <- (sim_dat$observed - min(sim_dat$observed))
  
  # Plot simulated "true" data
  print(ggplot(sim_dat, aes(X, Y)) +
    geom_tile(aes(fill = eta)) +
    facet_wrap(~year) +
    scale_color_gradient2()
  )
  #ggsave(paste0("true_density_range_",i,".png"))
  
  # Plot simulated "observed" data
  print(ggplot(sim_dat, aes(X, Y)) +
    geom_tile(aes(fill = observed)) +
    facet_wrap(~year) +
    scale_color_gradient2()
  )
 #ggsave(paste0("observed_density_range_",i,".png"))
  
#saveRDS(sim_dat, paste0("sim_dat_",i,".RDS"))
}


# Simulating from a fitted model ----

# fit model using example BC Pacific Cod bottom trawl data
mesh2 <- make_mesh(pcod, xy_cols = c("X", "Y"), cutoff = 10)

fit <- sdmTMB(
  density ~ s(depth),
  data = pcod,
  mesh = mesh2,
  family = tweedie(link = "log"),
  spatial = "on"
)

# predict to full domain
p_all <- predict(fit, newdata = qcs_grid)
ggplot(p_all, aes(X, Y, fill = exp(est))) + geom_raster() +
  scale_fill_viridis_c(trans = "sqrt")

# predict to fitted locations
p <- predict(fit)
ggplot(p, aes(X, Y, color = est)) + geom_point(aes(size = exp(est))) 

# simulate with default settings:
p$sim1 <- simulate(fit)
# simulate with new random fields:
p$sim2 <- simulate(fit, re_form = ~ 0)

ggplot(p, aes(X, Y, color = sim1)) + geom_point(aes(size = exp(sim1))) 
ggplot(p, aes(X, Y, color = sim2)) + geom_point(aes(size = exp(sim2))) 

# probably not the best example because fitted and simulated values differ a lot