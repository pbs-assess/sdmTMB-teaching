library(sdmTMB)
library(ggplot2)
options(ggplot2.continuous.colour = "viridis")
options(ggplot2.continuous.fill = "viridis")
theme_set(theme_minimal())

# ---------------------------------------------------------------------

# Try simulating some random fields to see how the range and sigma affect
# the random field.

# A grid of X and Y values:
predictor_dat <- expand.grid(
  X = seq(0, 1, length.out = 100), Y = seq(0, 1, length.out = 100)
)
mesh <- make_mesh(predictor_dat, c("X", "Y"), cutoff = 0.05)

# Re-run the following simulation and plotting repeatedly with
# different values for range and sigma_O. How do these affect
# the random field?
# Reminder: the range is the distance at which spatial correlation is
# effectively independent in units of X and Y.
# The X and Y values we're looking at scale from 0 to 1.

# run from here...
sim_dat <- sdmTMB_simulate(
  formula = ~ 1,
  data = predictor_dat,
  mesh = mesh,
  family = gaussian(link = "identity"),
  range = 0.6, # Try changing this (spatial range)
  sigma_O = 0.2, # Try changing this (spatial standard deviation)
  phi = 0.01, # observation error standard deviation; not used
  B = 0
)
ggplot(sim_dat, aes(X, Y, fill = mu)) +
  geom_raster() +
  scale_fill_gradient2()
# to here, repeatedly; try changing the range and sigma_O

# ---------------------------------------------------------------------

# Lets work with the built-in Pacific Cod data from British Columbia in
# Queen Charlotte Sound:
head(pcod)

# The density units are kg/km^2
# X and Y are coordinates in UTM zone 9
# We could have created those with `sdmTMB::add_utm_columns()`

# Look at the data:
ggplot(pcod, aes(X, Y, size = density)) + geom_point()

ggplot(pcod, aes(X, Y, size = density, colour = log(density + 1))) +
  geom_point() +
  facet_wrap(~year)

# Construct a mesh
# cutoff = 10 means the minimum length of triangle sides is 10 km:
mesh <- make_mesh(pcod, xy_cols = c("X", "Y"), cutoff = 10)
plot(mesh)

# ggplot alternative with inlabru:

# ggplot(pcod, aes(X, Y)) +
#   geom_point(aes(size = density)) +
#   inlabru::gg(mesh$mesh)

# Start by fitting a spatial model with a smoother for depth:

fit1 <- sdmTMB(
  density ~ s(depth),
  data = pcod,
  family = tweedie(link = "log"),
  mesh = mesh,
  spatial = "on",
  silent = FALSE
)

# Look at the maximum gradient and standard errors. Is this
# consistent with convergence?
max(fit1$gradients)
fit1$sd_report

# Inspect the model:
fit1

tidy(fit1, conf.int = TRUE)
tidy(fit1, effects = "ran_pars", conf.int = TRUE)

# Plot the conditional effect of depth by predicting across a gradient
# of depths. If we had other predictors in the model, we could set them
# to a fixed value such as their mean.
nd <- data.frame(depth = seq(0, 400, length.out = 100))

# Predict on our data frame:
p <- predict(
  fit1,
  newdata = nd,
  re_form = ~ 0, # means only include the fixed effects
  se_fit = TRUE # means calculate standard errors
)

# Plot it:
ggplot(p, aes(depth, exp(est),
  ymin = exp(est - 1.96 * est_se), ymax = exp(est + 1.96 * est_se))) +
  geom_line() + geom_ribbon(alpha = 0.4)

# Or try an experimental helper function for simple 1D smooths:
plot_smooth(fit1)
plot_smooth(fit1, ggplot = TRUE) + coord_cartesian(xlim = c(0, 400))

# Let's now predict on a grid that covers the entire survey.
# This dataframe `qcs_grid` is built into the package:
head(qcs_grid)
nd <- dplyr::filter(qcs_grid, year == 2003) # any one year will do
p <- predict(fit1, newdata = nd)

# Plot the two model components and their combined effects:

# Depth effect contribution:
# (Everything not a random field)
ggplot(p, aes(X, Y, fill = exp(est_non_rf))) +
  geom_raster() +
  coord_fixed()

# Spatial random field:
ggplot(p, aes(X, Y, fill = omega_s)) +
  geom_raster() +
  scale_fill_gradient2() +
  coord_fixed()

# Overall estimate in link (log) space:
ggplot(p, aes(X, Y, fill = est)) +
  geom_raster() +
  coord_fixed()

# Overall estimate:
ggplot(p, aes(X, Y, fill = exp(est))) +
  geom_raster() +
  coord_fixed()

# ---------------------------------------------------------------------

# Now, expand this model to include independent spatiotemporal random fields.
# To simplify this model, we will not include the effect of depth.
# The spatial random fields should largely account for this effect.
# We'll include mean estimates for each year with `0 + as.factor(year)`
# The `~ 0` part tells R to not include a separate intercept.
# Without it we'd be estimating the first year as an intercept and years
# after that as deviations from that.

fit2 <- sdmTMB(
  density ~ 0 + as.factor(year), #<<
  data = pcod,
  family = tweedie(link = "log"),
  mesh = mesh,
  spatial = "on",
  time = "year", #<<
  spatiotemporal = "iid", #<<
  silent = FALSE
)

fit2

# Inspect the sd_report and maximum gradient for fit2 -- has the model converged?

tidy(fit2, conf.int = TRUE)
tidy(fit2, effects = "ran_pars", conf.int = TRUE)

# What is larger: sigma_O (the spatial standard deviation) or 
# sigma_E (the spatiotemporal standard deviation)?

# Again, predict on the survey grid:
p2 <- predict(fit2, newdata = qcs_grid)

# Overall prediction in in log space:
ggplot(p2, aes(X, Y, fill = est)) +
  geom_raster() +
  facet_wrap(~year) +
  coord_fixed()

# What are these?
ggplot(p2, aes(X, Y, fill = exp(est_non_rf))) +
  geom_raster() +
  facet_wrap(~year) +
  coord_fixed()

# Those were the annual mean effects

# Why do these spatial random field effects look the same every year?
ggplot(p2, aes(X, Y, fill = omega_s)) +
  geom_raster() +
  facet_wrap(~year) +
  scale_fill_gradient2() +
  coord_fixed()

# What processes might this represent?

# Here are the spatiotemporal random field effects each year:
ggplot(p2, aes(X, Y, fill = epsilon_st)) +
  geom_raster() +
  facet_wrap(~year) +
  scale_fill_gradient2() +
  coord_fixed()

# What processes might this represent?

# ---------------------------------------------------------------------

# By default, sdmTMB shares the 'range' parameter across the spatial and
# spatiotemporal random fields because it's often hard to estimate.
#
# Reminder: the range is the distance at which spatial correlation is
# effectively independent (~ 0.13 correlation).
#
# Try fitting the same model but let the range be estimated
# separately for the spatial and spatiotemporal parts.
# Look at the help file: ?sdmTMB
# Hint, see the argument `share_range`

fit3 <- sdmTMB(
  density ~ 0 + as.factor(year),
  data = pcod,
  family = tweedie(link = "log"),
  mesh = mesh,
  spatial = "on",
  time = "year",
  spatiotemporal = "iid",
  share_range = FALSE, #<<
  silent = FALSE
)

fit3
tidy(fit3, "ran_pars", conf.int = TRUE)

# What is different now?
# Is the spatial or spatiotemporal range larger?
# What does this mean?

# ---------------------------------------------------------------------

# Let's calculate an area-weighted population index.

# All grid cells are 2 x 2 km = 4 km^2
pred2 <- predict(fit2, newdata = qcs_grid, return_tmb_object = TRUE, area = 4)
ind <- get_index(pred2, bias_correct = FALSE)

# In applied situations, you would want to set `bias_correct = TRUE`
# but this is much slower.
# Also see an experimental faster option ?get_index_sims

ggplot(ind, aes(year, est)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  ylim(0, NA)
