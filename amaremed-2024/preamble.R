xaringanthemer::style_mono_accent(
  base_color = "#202020",
  header_font_google = xaringanthemer::google_font("Raleway"),
  text_font_google = xaringanthemer::google_font("Open Sans"),
  code_font_google = xaringanthemer::google_font("Fira Mono"),
  title_slide_background_image = "images/logo-sdmTMB.png",
  title_slide_background_size = "14%",
  title_slide_background_position = "50% 90%",
  base_font_size = "20px",
  header_h1_font_size = "2.1rem",
  text_font_size = "1.5rem",
  code_font_size = "1.1rem",
  link_color = "#0047AB"
)

knitr_opts <- list(
  message = FALSE,
  warning = FALSE,
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  dpi = 300,
  out.width = "700px",
  fig.asp = 1 / 1.618,
  cache = TRUE,
  autodep = TRUE,
  cache.comments = TRUE,
  fig.align = "center",
  echo = FALSE
)

library(ggplot2)
library(dplyr)
library(sdmTMB)
options(ggplot2.continuous.colour = "viridis")
options(ggplot2.continuous.fill = "viridis")
ggplot2::theme_set(ggplot2::theme_minimal())
