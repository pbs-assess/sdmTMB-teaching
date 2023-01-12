d0 <- readRDS("~/src/gfsynopsis-2021/report/data-cache-april-2022/lingcod.rds")
d <- d0$survey_sets
d <- select(d, longitude, latitude, depth_m, catch_weight, doorspread_m,
  tow_length_m, duration_min, density_kgpm2, survey_abbrev, speed_mpm, year) |>
  filter(survey_abbrev == "SYN QCS")
sum(is.na(d$tow_length_m))
d <- mutate(d, area1 = doorspread_m * tow_length_m)
d <- mutate(d, area2 = doorspread_m * duration_min * speed_mpm)
d <- mutate(d, area_swept = ifelse(!is.na(area1), area1, area2))
d <- select(d, year, longitude, latitude, catch_weight, area_swept, density_kgpm2, depth_m) |>
  filter(!is.na(depth_m))
saveRDS(d, here::here("dfo-tesa-2023/data/lingcod-qcs.rds"))
