xaringanthemer::style_mono_accent(
  base_color = "#202020",
  header_font_google = xaringanthemer::google_font("Raleway"),
  text_font_google = xaringanthemer::google_font("Open Sans"),
  code_font_google = xaringanthemer::google_font("Fira Mono"),
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
