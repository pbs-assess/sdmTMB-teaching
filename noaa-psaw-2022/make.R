folder <- "noaa-psaw-2022"
files <- c(
  "01-intro-random-fields",
  "02-intro-sdmTMB",
  "03-tree-example",
  "04-pcod-example",
  "05-extra"
)

rm <- function(x) if (file.exists(x)) file.remove(x)
rm_folder <- function(x) if (file.exists(x)) unlink(x)
purrr::walk(files, function(.x) {
  cat(.x, "\n")
  f <- paste0(here::here(folder, .x), ".html")
  rm(f)
  f <- paste0(here::here(folder, .x), "_cache")
  rm_folder(f)
  f <- paste0(here::here(folder, .x), "_files")
  rm_folder(f)
})

# https://github.com/rstudio/rmarkdown/issues/1673
render_separately <- function(...) callr::r(
  function(...) rmarkdown::render(..., envir = globalenv()), args = list(...), show = TRUE)
purrr::walk(files, function(.x) {
  cat(.x, "\n")
  render_separately(paste0(here::here(folder, .x), ".Rmd"))
})
