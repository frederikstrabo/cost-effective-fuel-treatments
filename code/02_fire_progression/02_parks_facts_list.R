##########################################    02_parks_facts_list.R   ##########################################

################# Purpose: Create a list of all the MTBS fires in our sample 

################# Output: "FACTS_Parks_Fire_List.csv" saved in "data/intermediate"

rm(list=ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr,sf, tmap, magrittr, rnaturalearth, rnaturalearthdata, ggplot2, maps, lwgeom, rgeos, raster, stars, haven, stargazer, quantmod, lubridate, tidyr, ggpubr, 
               rgdal, exactextractr, tictoc, terra, gtools, here, fixest, modelsummary, readr, rdrobust, prism, parallel,tmaptools, 
               OpenStreetMap, maptiles, gifski)


# Set Path
here::i_am("code/MTBS_FACTS_Fires.R")

# Set CRS
CRS = 5070
# Load Functions 

source(here("code", "functions","tidy_facts.R"))
st_erase = function(x, y) st_difference(x, st_make_valid(st_union(st_combine(y)))) # Taken from here https://r-spatial.github.io/sf/reference/geos_binary_ops.html


########    Load in Datasets    ########

#### FACTS - Forest Service Fuel Treatments

facts <- st_read(here("data", "raw", "FACTS", "S_USA.Activity_HazFuelTrt_PL.shp"))
facts <- tidy_facts(facts, inf_yr = 2020, crs = CRS)
facts$ROW_ID <- 1:nrow(facts)

# distinguish completed & incomplete projects
facts_comp <- filter(facts, COMPLETE == 1)
facts_incomp <- filter(facts, COMPLETE == 0)

#### Load in Parks Daily Fire Perimeters

tif_files_1 <- list.files(path = here("data", "raw", "PARKS", "DOBs_from_Parks"), pattern = "\\.tif$", full.names = TRUE)
tif_files_2 <- list.files(path = here("data", "raw", "PARKS", "Jan2025_DOB_fires"), pattern = "\\.tif$", full.names = TRUE)
tif_files <- c(tif_files_1, tif_files_2)

# Create an empty list to store raster objects
raster_list <- list()
raster_stars_list <- list()

# Load each raster file into the list
for (file in tif_files) {
  raster_list[[file]] <- raster(file)
  raster_stars_list[[file]] <- read_stars(file)
}

#### Load in MTBS

# Burn Area Boundaries

mtbs <- st_read(here("data", "raw", "MTBS", "MTBS_Burn_Area", "S_USA.MTBS_BURN_AREA_BOUNDARY.shp")) %>%
  st_transform(mtbs, crs = 5070) %>%
  st_make_valid() %>%
  dplyr::rename(YEAR_MTBS = YEAR, MONTH_MTBS = STARTMONTH) 

mtbs <- filter(mtbs, FIRE_TYPE == "Unknown" | FIRE_TYPE == "Wildfire") # No Wildland Fire Use Fires


# Create a list of .tif fires that intersect with completed or planned fuel treatments

fires_int_df <- data.frame(matrix(ncol = 2, nrow = 0))

colnames(fires_int_df) <- c('ID', 'RAS_NUM')

for (i in 1:length(raster_list)){
  
  print(i)
  
  fire <- raster_list[[i]]
  
  # Get the source of the raster
  raster_source <- fire@file@name
  
  # Modify the source to remove "_1_dob.tif"
  fire_id <- sub("_1_dob\\.tif$", "", basename(raster_source))
  fire_id <- gsub("\\.dob.*$", "", fire_id)
  
  mtbs_fire <- filter(mtbs, FIRE_ID == fire_id)
  
  # intersect MTBS w/FACTS
  mtbs_facts_comp_int <- st_intersection(mtbs_fire, facts_comp)
  mtbs_facts_incomp_int <- st_intersection(mtbs_fire, facts_incomp)
  
  mtbs_facts_comp_int_filt <- as.data.frame(subset(mtbs_facts_comp_int, YEAR < YEAR_MTBS | YEAR == YEAR_MTBS & MONTH < MONTH_MTBS)) %>%
    filter(YEAR > 2004) %>%
    filter(YEAR_MTBS - YEAR < 10)
  
  mtbs_facts_incomp_int_filt <- as.data.frame(subset(mtbs_facts_incomp_int, YEAR_PLANNED < YEAR_MTBS | YEAR_PLANNED == YEAR_MTBS & MONTH_PLANNED < MONTH_MTBS)) %>%
    filter(YEAR_PLANNED > 2004)
  
  mtbs_facts_int_filt <- rbind(mtbs_facts_comp_int_filt, mtbs_facts_incomp_int_filt)
  
  nrow(mtbs_facts_int_filt) == 0
  
  if(nrow(mtbs_facts_int_filt) == 0){
    next
  }
  
  df_temp <- data.frame(matrix(ncol = 2, nrow = 0))
  
  colnames(df_temp) <- c('ID', 'RAS_NUM')
  df_temp[1, 1] <- fire_id
  df_temp[1, 2] <- i
  fires_int_df <- rbind(fires_int_df, df_temp)
  
}


## Save list of Parks fires that intersect as .csv

write_csv(fires_int_df, here("data", "intermediate", "FACTS_Parks_Fire_List.csv"))