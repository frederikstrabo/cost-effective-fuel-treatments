##########################################    01_define_fire_sample_mtbs_facts.R   ##########################################

################# Purpose: Intersect FACTS with MTBS to get a list of all the large wildfires that intersect with a completed or incomplete fuel treatments from USFS FACTS

################# Output: "FACTS_MTBS_Fire_List.csv" - saved in data\intermediate

################# Estimated run time: ~25 minutes 

rm(list=ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr,sf, tmap, magrittr, rnaturalearth, rnaturalearthdata, ggplot2, maps, lwgeom, rgeos, raster, stars, haven, stargazer, quantmod, lubridate, tidyr, ggpubr, 
               rgdal, exactextractr, tictoc, terra, gtools, here, fixest, modelsummary, readr, rdrobust, prism, parallel)


# Set Path
here::i_am("code/01_sample/01_define_fire_sample_mtbs_facts.R")
st_erase = function(x, y) st_difference(x, st_make_valid(st_union(st_combine(y)))) # Taken from here https://r-spatial.github.io/sf/reference/geos_binary_ops.html

# Load Functions 

source(here("code", "functions","tidy_facts.R"))

########    Load in Datasets    ########

#### FACTS - Forest Service Fuel Treatments

facts <- st_read(here("data", "raw", "FACTS", "S_USA.Activity_HazFuelTrt_PL.shp"))
facts <- tidy_facts(facts, inf_yr = 2023, crs = 5070)

# distinguish completed & incomplete projects
facts_comp <- filter(facts, COMPLETE == 1)
facts_incomp <- filter(facts, COMPLETE == 0)


# Load in MTBS

mtbs <- st_read(here("data", "raw", "MTBS", "MTBS_Burn_Area", "S_USA.MTBS_BURN_AREA_BOUNDARY.shp")) %>%
  st_transform(mtbs, crs = 5070) %>%
  st_make_valid()

mtbs <- filter(mtbs, FIRE_TYPE == "Unknown" | FIRE_TYPE == "Wildfire") # No Wildland Fire Use Fires

max(mtbs$YEAR)

# Load in USA

usa <- rnaturalearth::ne_states(country = "United States of America", returnclass = "sf") %>%
  st_transform(crs = 5070)

# 11 Western States - WA, OR, CA, NV, ID, MT, NM, AZ, UT, CO, WY

usa_west <- filter(usa, name == "Washington" | name == "Oregon" | name == "California" | name == "Nevada" | name == "Idaho" | 
                     name == "Montana" | name == "New Mexico" | name == "Arizona" | name == "Utah" | name == "Colorado" | name == "Wyoming") %>%
  st_transform(crs = 5070)


########    Define Sample of Fires  #######

## 1. Only look at MTBS fires that occur in the 11 Western U.S. States

mtbs_west <- st_crop(mtbs, st_bbox(usa_west))

## 2. Only look at fires that happened 2006-2023

mtbs_west <- filter(mtbs_west, YEAR > 2005 & YEAR <= 2023)

## 3. Only look at fires which intersect with a fuel treatment (complete or incomplete) that was complete or planned before the fire intersected with the treatment boundary

# Rename year & month to be different from facts
mtbs_west <- dplyr::rename(mtbs_west, YEAR_MTBS = YEAR, MONTH_MTBS = STARTMONTH) 

# intersect MTBS w/FACTS
mtbs_facts_comp_int <- st_intersection(mtbs_west, facts_comp)
mtbs_facts_incomp_int <- st_intersection(mtbs_west, facts_incomp)

mtbs_facts_comp_int_filt <- as.data.frame(subset(mtbs_facts_comp_int, YEAR < YEAR_MTBS | YEAR == YEAR_MTBS & MONTH < MONTH_MTBS)) %>%
  filter(YEAR > 2004) %>% 
  arrange(FIRE_ID)

mtbs_facts_incomp_int_filt <- as.data.frame(subset(mtbs_facts_incomp_int, YEAR_PLANNED < YEAR_MTBS | YEAR_PLANNED == YEAR_MTBS & MONTH_PLANNED < MONTH_MTBS)) %>%
  filter(YEAR_PLANNED > 2004) %>%
  arrange(FIRE_ID)

mtbs_facts_int <- rbind(mtbs_facts_comp_int_filt, mtbs_facts_incomp_int_filt)

MTBS_FACTS_List <- mtbs_facts_int %>%
  group_by(FIRE_ID) %>%
  summarise(FIRE_NAME = first(FIRE_NAME),
            IRWINID = first(IRWINID),
            YEAR = first(YEAR_MTBS),
            N_Treatments = n()) %>%
  arrange(desc(N_Treatments))

## Save list of fires as a .csv

write_csv(MTBS_FACTS_List, here("data", "intermediate", "FACTS_MTBS_Fire_List.csv"))