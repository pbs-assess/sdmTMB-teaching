PARALLEL <- TRUE

folder <- "noaa-psaw-2022"
files <- list.files(folder, pattern = "\\.Rmd$")
files <- gsub("\\.Rmd$", "", files)

rm <- function(x) if (file.exists(x)) file.remove(x)
rm_folder <- function(x) if (file.exists(x)) unlink(x, recursive = TRUE)
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
  function(...) rmarkdown::render(..., envir = globalenv()),
  args = list(...), show = TRUE)

if (!PARALLEL) {
  purrr::walk(files, function(.x) {
    cat(.x, "\n")
    render_separately(paste0(here::here(folder, .x), ".Rmd"))
  })
} else {
  future::plan(future::multisession)
  options(future.rng.onMisuse = "ignore")
  furrr::future_walk(files, function(.x) {
    render_separately(paste0(here::here(folder, .x), ".Rmd"))
  })
}

rmarkdown::render("noaa-psaw-2022/index.md")
