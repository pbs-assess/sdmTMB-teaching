theme_set(ggsidekick::theme_sleek())
library(ggplot2)

x <- seq(-1, 1, length.out = 100)
y <- seq(-1, 1, length.out = 100)
d <- expand.grid(x = x, y = y)
mesh <- make_mesh(d, xy_cols = c("x", "y"), cutoff = 0.05)

s <- sdmTMB_sim(
  x = d$x, y = d$y, mesh = mesh,
  time_steps = 6, rho = 0.8,
  phi = 0.2, range = 0.8, sigma_O = 0, sigma_E = 0.3,
  seed = 123, family = gaussian()
)
g <- ggplot(s, aes(x, y, fill = mu)) + geom_raster() +
  scale_fill_gradient2() +
  facet_wrap(~time) +
  coord_fixed(expand = FALSE)
ggsave("RF-rho-0.8.png", width = 8, height = 6)

s <- sdmTMB_sim(
  x = d$x, y = d$y, mesh = mesh,
  time_steps = 6, rho = 0,
  phi = 0.2, range = 0.8, sigma_O = 0, sigma_E = 0.3,
  seed = 123, family = gaussian()
)
g <- ggplot(s, aes(x, y, fill = mu)) + geom_raster() +
  scale_fill_gradient2() +
  facet_wrap(~time) +
  coord_fixed(expand = FALSE)
ggsave("RF-rho-0.png", width = 8, height = 6)
