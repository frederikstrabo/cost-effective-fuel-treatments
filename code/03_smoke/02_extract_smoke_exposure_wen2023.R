##########################################    02_extract_smoke_exposure_wen2023.R  ##########################################

################# Purpose: Uses data and code from Wen et al. (2023) to get population-day weighted PM2.5 smoke exposure estimates for fires in our sample 

#### Outputs: "Smoke_Fire_Effects_Wen2023.csv" saved in data/intermediate/Wen2023
####          - Each row is a unique MTBS fire with estimated deaths, death costs, earnings lost, and total population day weighted PM2.5 exposure

rm(list=ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr,sf, tmap, magrittr, rnaturalearth, rnaturalearthdata, ggplot2, maps, lwgeom, rgeos, raster, stars, haven, stargazer, quantmod, lubridate, tidyr, ggpubr, 
               rgdal, exactextractr, tictoc, terra, gtools, here, fixest, modelsummary, readr, rdrobust, prism, parallel,tmaptools, 
               OpenStreetMap, maptiles, gifski, purrr, nngeo, stringr, geosphere, data.table, fst, cowplot, pbmcapply, pbapply, patchwork, ggalluvial, ggpubr, RStoolbox)

# Set Path
here::i_am("code/03_smoke/02_extract_smoke_exposure_wen2023.R")


###### Load in Data sources

### population data ----
grid_pop_dt <- readRDS(here("data", "raw", "Wen2023", "clean", "population", "grid_pop_dt.rds"))
grid_pop_dt <- grid_pop_dt[, `:=`(year = as.numeric(as.character(year)))]

### globfire data ----
globfire_df <- st_read(here("data", "raw", "Wen2023", "clean", "globfire_na_final_area_2006-2020.shp"))

globfire_final_df <- globfire_df %>%
  arrange(IDate) %>%
  st_transform(crs = "epsg:5070") %>%
  mutate(
    fire_centroid_x = as.vector(sf::st_coordinates(st_centroid(st_geometry(.)))[, 1]), # store the centroid for dist calc later
    fire_centroid_y = as.vector(sf::st_coordinates(st_centroid(st_geometry(.)))[, 2])
  )


#### FRED Elderly Population Data

FRED_65_Over <- read_csv(here("data", "raw", "FRED", "SPPOP65UPTOZSUSA.csv")) %>% 
  dplyr::rename(Pct_Elderly_USA = SPPOP65UPTOZSUSA)

FRED_65_Over <- mutate(FRED_65_Over, YEAR = year(as.Date(observation_date))) %>% dplyr::select(YEAR, Pct_Elderly_USA)


#### MTBS Fires

mtbs <- st_read(here("data", "raw", "MTBS", "MTBS_Burn_Area", "S_USA.MTBS_BURN_AREA_BOUNDARY.shp")) %>%
  st_transform(mtbs, crs = 5070) %>%
  st_make_valid() %>%
  dplyr::rename(YEAR_MTBS = YEAR, MONTH_MTBS = STARTMONTH)

mtbs <- filter(mtbs, FIRE_TYPE == "Unknown" | FIRE_TYPE == "Wildfire") # No Wildland Fire Use Fires

fires_int_df <- read_csv(here("data", "intermediate", "FACTS_Parks_Fire_List.csv"))

MTBS_IDs <- unique(fires_int_df$ID)


#### Add extra fires 2007-2016 

fires_mtbs_int_df <- read_csv(here("data", "intermediate", "FACTS_MTBS_Fire_List.csv"))

# Filter to fires from 2007-2016

fires_mtbs_int_df_filt <- filter(fires_mtbs_int_df, YEAR %in% seq(2007, 2016, 1))

extra_fires <- unique(fires_mtbs_int_df_filt$FIRE_ID)

full_sample <- c(MTBS_IDs, extra_fires)

mtbs_sample <- filter(mtbs, FIRE_ID %in% full_sample)


#### Wen et al. 2023 code - goal is to get the population and day weighted PM2.5 smoke exposure 

## Figure 3 ----

## asthma ed estimates from heft-neal et al. 2023
ASTHMA_ED_COEFF <- 0.01712942 / 100000
ASTHMA_ED_SE <- 0.000828155 / 100000
ASTHMA_ED_BASE_RATE <- 0.9585968 / 100000

SAVE_PATH <- here("data", "raw", "Wen2023", "clean", "fire_smokepm")
WINDOW_SIZE <- 3

## iterate over year months to calculate total smoke at grid cells
MONTH_YEAR_LIST <- unique(zoo::as.yearmon(seq.Date(from = ymd("2006-04-01"), to = ymd("2020-12-31"), by = "days")))

total_smokepm_dt <- pbmclapply(MONTH_YEAR_LIST, function(MONTH_YEAR) {
  
  date_ <- as.Date(paste0(MONTH_YEAR, " 01"), format = "%b %Y %d")
  
  year_ <- year(date_)
  month_ <- str_pad(month(date_), 2, pad = "0")
  
  temp_fire_smokepm_dt <- read_fst(str_glue("{SAVE_PATH}/{WINDOW_SIZE}0km/{year_}-{month_}_{WINDOW_SIZE}0km_focal_grid_fire_smokepm_agg.fst"), as.data.table = T)
  
  ## aggregate smokePM to the fireID level (over days) then join pop data and aggregate to yearly smokepm over cells per fire
  temp_fire_smokepm_dt <- temp_fire_smokepm_dt[, `:=`(
    month = month(date),
    year = year(date)
  )][!is.na(contrib_smokePM),
     .(total_contrib_smokePM = sum(contrib_smokePM, na.rm = T)),
     by = c("ID", "fire_id", "month", "year")
  ][grid_pop_dt, `:=`(pop = i.pop), on = c("ID", "year")][, `:=`(
    pop_smokepm = pop * total_contrib_smokePM,
    asthma_ed_visits = pop * total_contrib_smokePM * ASTHMA_ED_COEFF,
    asthma_ed_visits_lb = pop * total_contrib_smokePM * (ASTHMA_ED_COEFF - qnorm(.95) * ASTHMA_ED_SE),
    asthma_ed_visits_ub = pop * total_contrib_smokePM * (ASTHMA_ED_COEFF + qnorm(.95) * ASTHMA_ED_SE)
  )][, .(
    affected_pop = sum(pop, na.rm = T),
    total_pop_smokePM = sum(pop_smokepm, na.rm = T),
    total_contrib_smokePM = sum(total_contrib_smokePM, na.rm = T),
    total_asthma_ed_visits = sum(asthma_ed_visits, na.rm = T),
    total_asthma_ed_visits_lb = sum(asthma_ed_visits_lb, na.rm = T),
    total_asthma_ed_visits_ub = sum(asthma_ed_visits_ub, na.rm = T)
  ),
  by = c("fire_id", "month", "year")
  ]
  
  return(temp_fire_smokepm_dt)
}) %>% bind_rows()

sorted_total_smokepm_dt <- total_smokepm_dt[, .(
  total_pop_smokePM = sum(total_pop_smokePM, na.rm = T),
  total_smokePM = sum(total_contrib_smokePM, na.rm = T),
  total_asthma_ed_visits = sum(total_asthma_ed_visits, na.rm = T),
  total_asthma_ed_visits_lb = sum(total_asthma_ed_visits_lb, na.rm = T),
  total_asthma_ed_visits_ub = sum(total_asthma_ed_visits_ub, na.rm = T)
),
by = c("fire_id")
][order(-total_pop_smokePM)]

total_smokepm_dt[, fire_id := as.numeric(fire_id)]

## Merge smoke exposure with Wen satelitte perimeters data

fire_smoke_merge <- merge(total_smokepm_dt, globfire_final_df, by.x = "fire_id", by.y = "Id", all.x = T)
fire_smoke_merge_filt <- filter(fire_smoke_merge, is.na(fire_id) == F)

fire_smoke_merge_filt_sf <- st_as_sf(fire_smoke_merge_filt, sf_column_name = "geometry") %>% st_transform(crs = 5070) %>% st_make_valid()

#### Intersect Wen satelitte perimeters with MTBS perimeters to know which fires are in our sample

fires_smoke_intersections <- st_intersection(mtbs_sample, fire_smoke_merge_filt_sf)[0,] %>% 
  dplyr::select(FIRE_ID, fire_id, affected_pop, total_pop_smokePM, total_contrib_smokePM, total_asthma_ed_visits)


for (fire in full_sample){
  
  print(fire)
  
  mtbs_fire <- filter(mtbs, FIRE_ID == fire)
  
  current_year <- mtbs_fire$YEAR_MTBS
  
  fire_smoke_merge_year <- filter(fire_smoke_merge_filt_sf, year == current_year)
  
  fires_smoke_int_temp <- st_intersection(mtbs_fire, fire_smoke_merge_year) %>% 
    dplyr::select(FIRE_ID, fire_id, affected_pop, total_pop_smokePM, total_contrib_smokePM, total_asthma_ed_visits)
  
  fires_smoke_intersections <- rbind(fires_smoke_intersections, fires_smoke_int_temp)
  
}


#### Calculate Fire Specific Effects

Fires_Smoke_Effects <- as.data.frame(fires_smoke_intersections) %>%
  group_by(FIRE_ID) %>%
  summarise(tot_affected_pop = sum(affected_pop), 
            total_pop_smokePM = sum(total_pop_smokePM), 
            total_asthma_ed_visits = sum(total_asthma_ed_visits))

mtbs_id_year <- dplyr::select(as.data.frame(mtbs_sample), FIRE_ID,  YEAR_MTBS)

Fires_Smoke_Effects <- merge(Fires_Smoke_Effects, mtbs_id_year, by = "FIRE_ID")

Fires_Smoke_Effects <- merge(Fires_Smoke_Effects, FRED_65_Over, by.x = "YEAR_MTBS", by.y = "YEAR")

#### Create estimated Death_f and EarningsLoss_f

## Deaths_f are calculated as smoke exposure weighted by per million people and 
#             the percent of the elderly multiplied by 0.69 estimate from Deryugina et al. 2019.
# DeathCosts are Deaths * VSL, which is 7.4 million USD (2006) taken from the EPA - https://www.epa.gov/environmental-economics/mortality-risk-valuation
#               we adjust 7.4 million to 11.16 million (2023 USD).

## EarningsLoss_f are calculated as smoke exposure weighted by the total days in a quarter - 92 (top end) 
#                     and multiplied by 103.1 from Borgshulte et al. 2022 we adjust this to 124.44 (2023 USD)

Fires_Smoke_Effects <- mutate(Fires_Smoke_Effects, 
                              Deaths = (total_pop_smokePM/1000000)*(Pct_Elderly_USA/100)*0.69,
                              EarningsLoss = (total_pop_smokePM/92)*124.44
)
Fires_Smoke_Effects <- mutate(Fires_Smoke_Effects, DeathCost = Deaths*11.16*1000000)



#### Save as .csv

Fires_Smoke_Effects <- dplyr::select(Fires_Smoke_Effects, FIRE_ID, Deaths, DeathCost, EarningsLoss, total_pop_smokePM)

dir.create(here("data", "intermediate", "Wen2023"))

write_csv(Fires_Smoke_Effects, here("data", "intermediate", "Wen2023", "Smoke_Fire_Effects_Wen2023.csv"))
