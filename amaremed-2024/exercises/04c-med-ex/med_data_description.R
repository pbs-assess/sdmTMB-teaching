# columns survey data:
#code: contains country, GSA number, year and haul, e.g: ITA_16_1999_1
#area: GSA number
#X: Longitude 
#Y: Latitude 
#year is the year of the survey
#month is the month of the survey in this year
#depth: depth in m
#swept.area: swept area in km^2 of the survey (in this year and month and haul)
#codeIndex: only one code used to merge survey data and grid data
#tmp_bot: bottom temperature in C°
#sal_bot: bottom salinity in psu
#tmp.sur: surface temperature in C° (first layer, ~ 1.43m)
#O2_bot: bottom oxygen in mmmo/m^3
#chl_int: Chlorophyll integrated 0-200m in mg/m^3
#poc_int: particulate organic carbon integrated 0-200m in mg/m^3
#poc_bot: particulate organic carbon at the bottom in mg/m^3
#NO3_bot: bottom nitrate mmol/m^3
#PO4_bot: phosphate at the bottom mmol/m^3
#SPECIE: Code to identify the species (genus and species, e.g. PAPE-LON)
#length_class: Length class of the individual measured in each haul
#nblon: Number of individuals measured for each length class, for each haul
#PA: Binary code for presence or absence, 0 or 1 respectively
#N_km2: Index of the number of individuals per km^2, calculated from nblon/sewpt.area for each haul
#ptot: kg of the species in the haul
#kg_km2: Index of kg/km^2, calculated from ptot/swetp.area (in this case the number is repeated for all nblon findings in the haul)
#X.utm: easting coordinates in UTM in m (zone 33)
#Y.utm: nothing coordinatwes in UTM in m (zone 33)

# columns grid:
#X: Longitude (resolution for grid 1/24°, ~ 4km)
#Y: Latitude (resolution for grid 1/24°, ~ 4km)3
#years of the time series
#month we used the month of July, the most common month for survey data (in this area, see row 55-56)
#index: only one code used to merge survey data and grid data
#depth: depth in m of the Copernicus grid (1/24°)
#sal_bot: bottom salinity in psu
#tmp.sur: surface temperature in C° (first layer, ~ 1.43m)
#tmp_bot: bottom temperature in C°
#O2_bot: bottom oxygen in mmmo/m^3
#chl_int: Chlorophyll integrated 0-200m in mg/m^3
#poc_int: particulate organic carbon integrated 0-200m in mg/m^3
#poc_bot: particulate organic carbon at the bottom in mg/m^3
#NO3_bot: bottom nitrate mmol/m^3
#PO4_bot: phosphate at the bottom mmol/m^3
#X.utm: easting coordinates in UTM in m (zone 33)
#Y.utm: nothing coordinates in UTM in m (zone 33)


# 'survey data'
# merl_TATBTC.AMMED <- readRDS('merl_TATBTC.AMMED.rds')
# pape_TATBTC.AMMED <- readRDS('pape_TATBTC.AMMED.rds')
# 
# 'grid'
# table(merl_TATBTC.AMMED$month)
# table(pape_TATBTC.AMMED$month)
# 
# grid_1516 <- readRDS('grid_1516.rds')





