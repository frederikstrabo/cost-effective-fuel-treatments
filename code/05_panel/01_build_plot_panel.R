##########################################    01_build_plot_panel.R  ##########################################

################# Purpose: This code develops the function "create_SpatialDiD_plot_panel" which for each fire 
#################           in our sample creates plots, defined by their unique direction and distance from ignition.
#################           For each plot we calculate their associated wildfire outcomes, treatment status, and relevant 
#################           covariates to be used in our spatial DiD analysis. The function allows for flexibility to change 
#################           the number of radial directions = L, the distance of each plot = K, both of which will be tested in "07_robustness".

################  ***Note***: This code defaults to running the baseline sample. At the end of the script the user can uncomment code used to construct
###############               alternative samples used in "07_robustness", which greatly increases the run time.

################# Output: "SpatialDiD_Grids_L24_K05.csv" in the folder "data/temp" - the panel is finished in "06_analysis/01_descriptive_stats.R" where smoke exposure is imputed for fires without PM2.5 exposure estimates

################# Optional outputs: i) Different direction panels: "SpatialDiD_Grids_L36_K05.csv", "SpatialDiD_Grids_L18_K05.csv", and "SpatialDiD_Grids_L12_K05.csv". 
#################                   ii) reported ignition points: "SpatialDiD_Grids_L24_K05_V2.csv".
#################                   iii) alternative treatment thresholds: "SpatialDiD_Grids_L24_K05_100p.csv", "SpatialDiD_Grids_L24_K05_75p.csv", "SpatialDiD_Grids_L24_K05_25p.csv", "SpatialDiD_Grids_L24_K05_zero.csv".
#################                   all saved in "data/intermediate"


################# Estimated run time: ~ 4 hours & 15 mins for just baseline panel & ~34 hours for the baseline panel & alternative samples

rm(list=ls())

# if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr,sf, tmap, magrittr, rnaturalearth, rnaturalearthdata, ggplot2, maps, lwgeom, raster, stars, haven, stargazer, quantmod, lubridate, tidyr, ggpubr, 
               exactextractr, tictoc, terra, gtools, here, fixest, modelsummary, readr, rdrobust, prism, parallel,tmaptools, 
               OpenStreetMap, maptiles, gifski, purrr, nngeo, stringr, geosphere)

tic()

# Set Path
here::i_am("code/05_panel/01_build_plot_panel.R")

# Load Functions 

source(here("code", "functions","tidy_facts.R"))
source(here("code", "functions","calculate_distance.R"))
st_erase = function(x, y) st_difference(x, st_make_valid(st_union(st_combine(y)))) # Taken from here https://r-spatial.github.io/sf/reference/geos_binary_ops.html

########    Load in Datasets    ########

#### FACTS - Forest Service Fuel Treatments

facts <- st_read(here("data", "raw", "FACTS", "S_USA.Activity_HazFuelTrt_PL.shp"))
facts <- tidy_facts(facts, inf_yr = 2023, crs = 5070)
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

fire_lookup <- tibble(
  file = tif_files,
  filename = basename(tif_files),
  ID = filename |>
    str_remove("\\.tif$") |>
    str_remove("\\.dob$") |>
    str_remove("_1_dob$") |>
    str_remove("_dob$")
)

#### Load in USFS NF & Wilderness Areas

USFS_NF <- st_read(here("data", "raw", "USFS", "NationalForests", "S_USA.AdministrativeForest.shp")) %>% st_transform(crs = 5070)
WAs <- st_read(here("data", "raw", "USFS", "WildernessAreas", "S_USA.Wilderness.shp"))  %>% st_transform(crs = 5070)


#### Load in MTBS

# Burn Area Boundaries

mtbs <- st_read(here("data", "raw", "MTBS", "MTBS_Burn_Area", "S_USA.MTBS_BURN_AREA_BOUNDARY.shp")) %>%
  st_transform(mtbs, crs = 5070) %>%
  st_make_valid() %>%
  dplyr::rename(YEAR_MTBS = YEAR, MONTH_MTBS = STARTMONTH)

mtbs <- filter(mtbs, FIRE_TYPE == "Unknown" | FIRE_TYPE == "Wildfire") # No Wildland Fire Use Fires

mtbs_point <- st_read(here("data", "raw", "MTBS", "mtbs_fod_pts_data", "mtbs_FODpoints_DD.shp")) %>%
  st_transform(mtbs, crs = 5070) %>%
  st_make_valid()  %>%
  dplyr::rename(FIRE_ID = Event_ID)

#### Fire Suppression Lines

QAQC <- st_read(here("data", "raw", "QAQC", "FTE_containment.shp")) %>% st_transform(crs = 5070)
qaqc_ids <- unique(QAQC$MTBS_id)
mtbs_qaqc <- filter(mtbs, FIRE_ID %in% qaqc_ids)

#### Load in list of fires that intersect with FTs

fires_int_df <- read_csv(here("data", "intermediate", "FACTS_Parks_Fire_List.csv"))

#### Load WFEIS Fire Smoke Emissions

smoke_data <- read_csv(here("data", "raw", "WFEIS", "WFEIS_data.csv"))

#### Load Wen et al. (2023) Fire Specific PM2.5 effects

smoke_costs <- read_csv(here("data", "intermediate", "Wen2023", "Smoke_Fire_Effects_Wen2023.csv"))

#### Load in USFS Roads & US Highways

USFS_Roads <- st_read(here("data", "raw", "USFS", "Roads", "S_USA.RoadCore_FS.shp")) %>% st_transform(crs = 5070)
US_Highways <- st_read(here("data", "raw", "Highways", "tl_2016_us_primaryroads.shp")) %>% st_transform(crs = 5070)
Western_WUI <- st_read(here("data", "intermediate", "WUI", "western_wui_2010.shp")) %>% st_transform(crs = 5070)

###### Fuel Type & Fire Regime Characteristics (MFRI)

MFRI <- raster(here("data", "raw", "LANDFIRE", "US_140_MFRI", "Tif", "us_140mfri.tif"))
FBFM13 <- raster(here("data", "raw", "LANDFIRE", "US_105_FBFM13", "Tif", "us_105fbfm13.tif")) # Fuel Model Type in LANDFIRE 2001


###### Load topographic LANDFIRE data

# Load Slope percent, aspect, & Elev

slope <- raster(here("data", "raw", "LANDFIRE", "LF2020_SlpD_220_CONUS", "Tif", "LC20_SlpD_220.tif"))
aspect <- raster(here("data", "raw", "LANDFIRE", "LF2020_Asp_220_CONUS", "Tif", "LC20_Asp_220.tif"))
elev <- raster(here("data", "raw", "LANDFIRE", "LF2020_Elev_220_CONUS", "Tif", "LC20_Elev_220.tif"))
TRI <- terrain(elev, opt = "TRI")
top_landfire <- stack(slope, aspect, elev, TRI) %>% terra::rast()

## Houses & Structures

Buildings <- raster(here("data", "raw", "CommunitiesRisk", "CONUS", "BuildingCount_CONUS.tif"))
HUCount <- raster(here("data", "raw", "CommunitiesRisk", "CONUS", "HUCount_CONUS.tif"))


###### LAT Drops 

LAT_2017 <- st_read(here("data", "raw", "LAT", "drops17.shp")) %>% st_transform(crs = 5070)
LAT_2017$YEAR <- 2017
LAT_2017 <- dplyr::select(LAT_2017, YEAR)

LAT_2018 <- st_read(here("data", "raw", "LAT", "drops18.shp")) %>% st_transform(crs = 5070)
LAT_2018$YEAR <- 2018
LAT_2018 <- dplyr::select(LAT_2018, YEAR)

LAT_2019 <- st_read(here("data", "raw", "LAT", "drops19.shp")) %>% st_transform(crs = 5070)
LAT_2019$YEAR <- 2019
LAT_2019 <- dplyr::select(LAT_2019, YEAR)

LAT_2020_2021 <- st_read(here("data", "raw", "LAT", "drops20_21.shp")) %>% st_transform(crs = 5070)
LAT_2020_2021 <- mutate(LAT_2020_2021, YEAR = localYear1)
LAT_2020_2021 <- dplyr::select(LAT_2020_2021, YEAR)

LAT_2022 <- st_read(here("data", "raw", "LAT", "2022_Full_Year_03112024_VLAT.shp")) %>% st_transform(crs = 5070)
LAT_2022 <- mutate(LAT_2022, YEAR = Year)
LAT_2022 <- dplyr::select(LAT_2022, YEAR)

LAT_2023 <- st_read(here("data", "raw", "LAT", "2023_Full_Year_03112024_VLAT.shp")) %>% st_transform(crs = 5070)
LAT_2023 <- mutate(LAT_2023, YEAR = Year)
LAT_2023 <- dplyr::select(LAT_2023, YEAR)

LAT <- rbind(LAT_2017, LAT_2018, LAT_2019, LAT_2020_2021, LAT_2022, LAT_2023)
  
#### Load ACS Data

## Shapefiles of ACS Census blocks

ACS_sf <- st_read(here("data", "raw", "ACS", "nhgis0003_shape", "US_blck_grp_2020.shp"))
ACS_sf <- dplyr::select(ACS_sf, GISJOIN) %>% st_transform(crs = 5070)

## Create ACS sfs for 2013-2017, ..., 2018-2022 Community Surveys

ACS_2018_2022 <- read_csv(here("data", "raw", "ACS", "nhgis0003_csv", "nhgis0003_ds262_20225_blck_grp.csv"))
ACS_2018_2022 <- dplyr::select(ACS_2018_2022, GISJOIN, AQU4E001) %>% dplyr::rename(MED_HOUSE_VAL_2022 = AQU4E001)

ACS_2017_2021 <- read_csv(here("data", "raw", "ACS", "nhgis0003_csv", "nhgis0003_ds254_20215_blck_grp.csv"))
ACS_2017_2021 <- dplyr::select(ACS_2017_2021, GISJOIN, AOULE001) %>% dplyr::rename(MED_HOUSE_VAL_2021 = AOULE001)

ACS_2016_2020 <- read_csv(here("data", "raw", "ACS", "nhgis0003_csv", "nhgis0003_ds249_20205_blck_grp.csv"))
ACS_2016_2020 <- dplyr::select(ACS_2016_2020, GISJOIN, AMWBE001) %>% dplyr::rename(MED_HOUSE_VAL_2020 = AMWBE001)

ACS_2015_2019 <- read_csv(here("data", "raw", "ACS", "nhgis0003_csv", "nhgis0003_ds244_20195_blck_grp.csv"))
ACS_2015_2019 <- dplyr::select(ACS_2015_2019, GISJOIN, AL1HE001) %>% dplyr::rename(MED_HOUSE_VAL_2019 = AL1HE001) %>% 
  mutate(MED_HOUSE_VAL_2019 = as.numeric(MED_HOUSE_VAL_2019))

ACS_2014_2018 <- read_csv(here("data", "raw", "ACS", "nhgis0003_csv", "nhgis0003_ds239_20185_blck_grp.csv"))
ACS_2014_2018 <- dplyr::select(ACS_2014_2018, GISJOIN, AJ3QE001) %>% dplyr::rename(MED_HOUSE_VAL_2018 = AJ3QE001) %>% 
  mutate(MED_HOUSE_VAL_2018 = as.numeric(MED_HOUSE_VAL_2018))

ACS_2013_2017 <- read_csv(here("data", "raw", "ACS", "nhgis0003_csv", "nhgis0003_ds233_20175_blck_grp.csv"))
ACS_2013_2017 <- dplyr::select(ACS_2013_2017, GISJOIN, AH53E001) %>% dplyr::rename(MED_HOUSE_VAL_2017 = AH53E001) %>% 
  mutate(MED_HOUSE_VAL_2017 = as.numeric(MED_HOUSE_VAL_2017))

ACS_sf <- merge(ACS_sf, ACS_2018_2022, by = "GISJOIN", all.x = T) 
ACS_sf <- merge(ACS_sf, ACS_2017_2021, by = "GISJOIN", all.x = T)
ACS_sf <- merge(ACS_sf, ACS_2016_2020, by = "GISJOIN", all.x = T) 
ACS_sf <- merge(ACS_sf, ACS_2015_2019, by = "GISJOIN", all.x = T) 
ACS_sf <- merge(ACS_sf, ACS_2014_2018, by = "GISJOIN", all.x = T) 
ACS_sf <- merge(ACS_sf, ACS_2013_2017, by = "GISJOIN", all.x = T)

# Inflation adjustment
cpi_data <- read_csv(here("data", "raw", "FRED", "CPIAUCSL.csv")) %>%
  mutate(YEAR = year(observation_date))

# Average CPI by year
avg_cpi <- cpi_data %>%
  group_by(YEAR) %>%
  summarise(AVG_CPI = mean(CPIAUCSL, na.rm = TRUE))

# Reference CPI (2023)
ref_cpi <- avg_cpi %>%
  filter(YEAR == 2023) %>%
  pull(as.numeric(AVG_CPI))

# Inflation adjustment factor
avg_cpi <- avg_cpi %>%
  mutate(CPI_FACTOR = AVG_CPI / ref_cpi) %>%
  dplyr::select(YEAR, CPI_FACTOR)

adj_2023 <- avg_cpi[avg_cpi$YEAR == 2023, ]$CPI_FACTOR
adj_2022 <- avg_cpi[avg_cpi$YEAR == 2022, ]$CPI_FACTOR
adj_2021 <- avg_cpi[avg_cpi$YEAR == 2021, ]$CPI_FACTOR
adj_2020 <- avg_cpi[avg_cpi$YEAR == 2020, ]$CPI_FACTOR
adj_2019 <- avg_cpi[avg_cpi$YEAR == 2019, ]$CPI_FACTOR
adj_2018 <- avg_cpi[avg_cpi$YEAR == 2018, ]$CPI_FACTOR
adj_2017 <- avg_cpi[avg_cpi$YEAR == 2017, ]$CPI_FACTOR

## Adjust median housing values for inflation

ACS_sf <- mutate(ACS_sf, MED_HOUSE_VAL_2022 = MED_HOUSE_VAL_2022/adj_2022,
                 MED_HOUSE_VAL_2021 = MED_HOUSE_VAL_2021/adj_2021,
                 MED_HOUSE_VAL_2020 = MED_HOUSE_VAL_2020/adj_2020,
                 MED_HOUSE_VAL_2019 = MED_HOUSE_VAL_2019/adj_2019,
                 MED_HOUSE_VAL_2018 = MED_HOUSE_VAL_2018/adj_2018,
                 MED_HOUSE_VAL_2017 = MED_HOUSE_VAL_2017/adj_2017)


###### Step 1. Construct set of grids to use for Spatial RD Analysis

#### The function "create_SpatialDiD_plot_panel": 
####                                        1) Constructs plots for each fire in our sample defined by their direction and distance from the ignition point.
####                                        2) 

#### Relevant parameters: L = number of directions, K = length of distance bin in km, Tau = percent of a plot that needs to be treated to be considered "Treated", 
####    CENT is an indicator equal to 1 if ignition points are imputed by the centroid of the first day of burning polygon otherwise it defaults to using reported ignition points.

## recommended to first check if the code works by manually running the script on a single fire i.e. i = 1

L = 24
K = 0.5
Tau = 0.5
CENT = T


create_SpatialDiD_plot_panel <- function(L, K, Tau, CENT){
  
  #### Initialize our panel
  
  panel <- data.frame(matrix(ncol = 58, nrow = 0), geometry = st_sfc())
  
  colnames(panel) <- c("ID", "direction", "distance_bin", "TREAT_ID", "bearing", "USFS_NF", "Wilderness", "Grid_Acres", "Pct_Treated", "Treated", "COMPLETE",
                                  "TREAT_COST","Treated_LAG", "L_TAU", "Treated_Dir", "AFT_TREAT", "PRE", "POST", "time_to_treat", "TREAT_CAT", "TREAT_SIZE", "TREAT_TYPES", "DAY_BURN", 
                                  "PCT_BURN", "BURN", "BURN_LAG", "YEAR", "FIRE_SIZE", "YEARS_SINCE_TREAT", "DAY_BURN_CLOSEST", "DATE_BURNED", "BURN_SEV", "BURN_SEV0", 
                                  "BURN_SEV_MED_HIGH_PCT","BURN_SEV_HIGH_PCT",
                                  "FIRE_ID", "FIRE_NAME", "HIT_TREAT","DIST_TREAT", "DIST_HIT_TREATS", "Pct_Burned_10Y", "PREV_BURN_10Y", 
                                  "direction_distance_fire_FE", "direction_fire_FE", "distance_fire_FE", "Burn_Right_Away", "MTT", "MTT_LAG", "DELTA_MTT", "NA_DELTA_MTT",
                                  "FIRE_INTENSITY", "LOG_FIRE_INTENSITY", "SUP_LINES_PRESENT", "SUP_LINE_DIST", "SUP_LINE", "LOG_SUP_LINE_DIST", "SUP_LINE_LEN", "SUP_LINE_INT",
                                  "geometry")
  
  
  panel <- st_as_sf(panel, crs = 5070)
  
  #### Step 1. Loop through each fire in our sample, construct the plots, relevant wildfire outcomes, treatment status and other covariates.
  
  tic()
  
  ## Loop takes about 3 hours & 20 mins to run
  
  # L is the number of direction, K is the distance in kilometers of each distance bin
  for (i in 1:nrow(fires_int_df)){
    
    print(i)
    
    n = i
    
    ###### Step 1. get relevant fire characteristics, ignition point, load in relevant fire specific data, and intersect with treatments
    
    # Grab fire and raster ID
    fire_id <- as.character(fires_int_df[i, 1])
    ras_id <- as.numeric(fires_int_df[i, 2])
    
    # Get DOB perimeter, MTBS perimeter, fire name, year and size
    fire <- raster_list[[ras_id]]
    mtbs_fire <- filter(mtbs, FIRE_ID == fire_id)
    fire_name <- mtbs_fire$FIRE_NAME
    fire_year <- mtbs_fire$YEAR_MTBS
    fire_size <- mtbs_fire$ACRES
    
    #### Make sure the day of burning raster & MTBS fire ID match up 
    
    ras_fire_id <- fire_lookup[ras_id , "ID"]
    
    if (ras_fire_id != fire_id) {
      stop(
        paste0(
          "ERROR: MTBS fire ID mismatch.\n",
          "  Expected MTBS ID: ", mtbs_id, "\n",
          "  Raster layer name: ", raster_id, "\n",
          "Please check the file ordering or naming of the day-of-burning rasters in data/raw/PARKS."
        ),
        call. = FALSE
      )
    }
    
    # reported ignition point
    mtbs_fire_ig <- filter(mtbs_point, FIRE_ID == fire_id)
    
    # load MTT outputs
    MTT <- raster(here("data", "raw", "FB", "TestMTT", "MTT_Inputs", 
                       "MTT_Output", paste0("fire", n, "_ArrivalTime.asc")))
    
    MTT_INT <- raster(here("data", "raw", "FB", "TestMTT", "MTT_Inputs", 
                           "MTT_Output", paste0("fire", n, "_INTENSITY.asc")))
    
    #### Get ignition point from centroid of first day
    
    if (CENT == T){
      
      # Step 1: Find the minimum value in the raster
      min_day <- minValue(fire)
      
      # Step 2: Create a mask for cells with the minimum value
      ignition_mask <- fire == min_day
      
      # Step 3: Convert the mask to polygons
      ignition_polygons <- rasterToPolygons(ignition_mask, fun = function(x) x == 1, dissolve = TRUE)
      
      # Step 4: Convert polygons to an sf object for centroids
      ignition_sf <- st_as_sf(ignition_polygons)
      
      # Step 5: Compute the centroid
      mtbs_fire_ig <- st_centroid(ignition_sf) %>% st_transform(crs = 5070)
      
    }
    
    # load in MTBS burn severity raster
    mtbsBS_year <- raster(here("data", "raw", "MTBS", "MTBS_BSmosaics", fire_year, paste0("mtbs_CONUS_", fire_year, ".tif")))
    
    # intersect MTBS w/FACTS - both complete and incomplete treatments
    mtbs_facts_comp_int <- st_intersection(mtbs_fire, facts_comp)
    mtbs_facts_incomp_int <- st_intersection(mtbs_fire, facts_incomp)
    
    # filter to only treatments that were complete before the fire
    mtbs_facts_comp_int_filt <- as.data.frame(subset(mtbs_facts_comp_int, YEAR < YEAR_MTBS | YEAR == YEAR_MTBS & MONTH < MONTH_MTBS)) %>%
      filter(YEAR > 2004) %>%
      filter(YEAR_MTBS - YEAR <= 10)
    
    mtbs_facts_incomp_int_filt <- as.data.frame(subset(mtbs_facts_incomp_int, YEAR_PLANNED < YEAR_MTBS | YEAR_PLANNED == YEAR_MTBS & MONTH_PLANNED < MONTH_MTBS)) %>%
      filter(YEAR_PLANNED > 2004)
    
    mtbs_facts_int_filt <- rbind(mtbs_facts_comp_int_filt, mtbs_facts_incomp_int_filt)
    
    # grab treatment IDs
    treatment_IDs <- unique(mtbs_facts_int_filt$ROW_ID)
    
    # the bounding perimeter of the fire
    fire_extent_poly <- st_as_sfc(st_bbox(fire))
    
    # Find intersecting fuel treatments - subset to only complete treatments
    intersecting_facts <- filter(facts, ROW_ID %in% treatment_IDs)
    intersecting_facts <- filter(intersecting_facts, COMPLETE == 1) # Only look at complete projects for now
    
    # if no treatments skip this round in the loop - occurs when only looking at complete treatments
    if (nrow(intersecting_facts) == 0){
      next
    }
    
    ###### Step 2. Construct plots by their unique direction and distance from the ignition point
    
    # Get the boundary of the fire perimeter
    fire_boundary <- st_boundary(mtbs_fire)
    
    # Convert fire boundary to multipoint for easier distance calculations
    boundary_points <- st_cast(fire_boundary, "POINT")
    
    # Calculate the distances from the ignition point to each boundary point
    distances <- st_distance(mtbs_fire_ig, boundary_points)
    
    # Find the maximum distance - this sets the maximum distance bin used to construct our plots
    max_distance <- max(as.numeric(distances))  # Extract numeric value from the distance object
    
    # add 3 kilometers to max distance - allows for simulating fire growth in the absence of treatment
    max_distance <- max_distance + 3000
    
    # Define the number of directions
    num_rays <- L
    
    # Create a sequence of angles (in radians) that defines the directions of the fire
    angles <- seq(0, 2 * pi, length.out = num_rays + 1)[-1]
    
    # Function to calculate the end point of a ray given an origin, angle, and length
    calc_end_point <- function(origin, angle, length) {
      coords <- st_coordinates(origin)
      x <- coords[1]
      y <- coords[2]
      x_end <- x + length * cos(angle)
      y_end <- y + length * sin(angle)
      st_sfc(st_point(c(x_end, y_end)), crs = st_crs(origin))
    }
    
    # Calculate end points using the maximum distance
    end_points_max <- lapply(angles, calc_end_point, origin = mtbs_fire_ig, length = max_distance) %>%
      do.call(c, .)
    
    # Create rays using the maximum distance
    final_rays <- st_sfc(lapply(end_points_max, function(end_point) {
      st_linestring(rbind(st_coordinates(mtbs_fire_ig), st_coordinates(end_point)))
    }), crs = st_crs(mtbs_fire_ig))
    
    # Combine into a single sf object
    rays_sf <- st_sf(geometry = final_rays)
    rays_buffer <- st_buffer(mtbs_fire_ig, dist = max_distance)
    
    # Define the number of angles (e.g., 24 by default, but this can be changed)
    num_angles <- num_rays
    
    # Create a sequence of angles (in radians) based on the number of angles
    angles <- seq(0, 2 * pi, length.out = num_angles + 1)[-1]
    
    # Define distance bins (0.5 km intervals)
    distance_bins <- seq(0, max_distance, by = K*1000)  # K = 0.5  then 0.5*1000 meters = 0.5 km
    
    # Function to create wedges (pie slices)
    create_wedge_polygon <- function(angle1, angle2, distance_start, distance_end, origin) {
      # Convert angles to x, y coordinates for the end points of the rays
      coords <- st_coordinates(origin)
      x <- coords[1]
      y <- coords[2]
      
      # Ray 1 end point
      x1_start <- x + distance_start * cos(angle1)
      y1_start <- y + distance_start * sin(angle1)
      x1_end <- x + distance_end * cos(angle1)
      y1_end <- y + distance_end * sin(angle1)
      
      # Ray 2 end point
      x2_start <- x + distance_start * cos(angle2)
      y2_start <- y + distance_start * sin(angle2)
      x2_end <- x + distance_end * cos(angle2)
      y2_end <- y + distance_end * sin(angle2)
      
      # Create the polygon for this wedge
      polygon_coords <- rbind(
        c(x1_start, y1_start),  # start of ray 1
        c(x1_end, y1_end),      # end of ray 1
        c(x2_end, y2_end),      # end of ray 2
        c(x2_start, y2_start),  # start of ray 2
        c(x1_start, y1_start)   # close the polygon
      )
      
      st_sfc(st_polygon(list(polygon_coords)), crs = st_crs(origin))
    }
    
    # Initialize a list to store all direction-distance polygons
    direction_distance_polygons <- list()
    
    # Iterate through each angle (based on num_angles) and each distance bin
    for (i in 1:num_angles) {
      angle1 <- angles[i]
      angle2 <- if (i == num_angles) angles[1] else angles[i + 1]  # Wrap around at the last angle
      
      for (j in seq_along(distance_bins)[-1]) {  # Skip the first distance bin (0)
        # Get the start and end of the current distance bin
        bin_start <- distance_bins[j - 1]
        bin_end <- distance_bins[j]
        
        # Create the wedge polygon for the current direction-distance bin
        wedge_polygon <- create_wedge_polygon(angle1, angle2, bin_start, bin_end, mtbs_fire_ig)
        
        # Add metadata (direction and distance bin)
        direction_distance_polygons[[length(direction_distance_polygons) + 1]] <- st_sf(
          geometry = wedge_polygon,
          direction = i,  # Ray index as direction
          distance_bin_start = bin_start,
          distance_bin_end = bin_end
        )
      }
    }
    
    # Combine all polygons into one sf object
    direction_distance_sf <- do.call(rbind, direction_distance_polygons)
    direction_distance_sf <- st_make_valid(direction_distance_sf)
    
    
    #### Get the dominant angle of spread for each direction - used when creating control variable of wind difference
    
    # Ensure fire origin is a POINT and extract its coordinates
    fire_origin_coords <- st_coordinates(st_transform(mtbs_fire_ig, crs = 4326))  # (X = longitude, Y = latitude)
    
    direction_distance_sf <- st_transform(direction_distance_sf, crs = 4326)
    
    # Compute bearing for each polygon slice
    direction_distance_sf$bearing <- sapply(st_geometry(direction_distance_sf), function(geom) {
      # Get centroid of the polygon slice
      slice_centroid <- st_centroid(geom)
      slice_coords <- st_coordinates(slice_centroid)  # Extract coords
      
      # Compute the bearing from fire origin to slice centroid
      bearing_value <- bearing(c(fire_origin_coords[1], fire_origin_coords[2]), 
                               c(slice_coords[1], slice_coords[2]))
      
      # Ensure bearing is within 0-360°
      if (bearing_value < 0) bearing_value <- bearing_value + 360
      
      return(bearing_value)
    })
    
    direction_distance_sf <- st_transform(direction_distance_sf, crs = 5070)
    
    ###### Step 3. With our now defined plots extract relevant controls, wildfire outcomes and treatment status
    
    plots <- direction_distance_sf
    plots$Plot_ID <- 1:nrow(plots)
    plots_Cent <- st_centroid(plots)
    
    # Calculate acres in a plot
    plots$Grid_Acres <- as.numeric(st_area(plots)/4046.86)
    
    ###### Create National Forest & Wilderness indicators
    
    ## Intersect with National Forests
    plots_NF_int <- st_intersection(plots_Cent, USFS_NF)
    
    plots_usfs_nf_int <- as.data.frame(plots_NF_int) %>%
      group_by(Plot_ID) %>%
      summarise(USFS_NF = 1) %>%
      dplyr::select(Plot_ID, USFS_NF)
    
    plots <- merge(plots, plots_usfs_nf_int, by = "Plot_ID", all.x = T)
    plots[is.na(plots$USFS_NF), "USFS_NF"] <- 0
    
    ## Intersect with Wilderness Areas
    plots_NF_int <- st_intersection(plots_Cent, WAs)
    
    plots_usfs_wa_int <- as.data.frame(plots_NF_int) %>%
      group_by(Plot_ID) %>%
      summarise(Wilderness = 1) %>%
      dplyr::select(Plot_ID, Wilderness)
    
    plots <- merge(plots, plots_usfs_wa_int, by = "Plot_ID", all.x = T)
    plots[is.na(plots$Wilderness), "Wilderness"] <- 0
    
    
    ###### Extract relevant treatment information for each plot 
    
    # Group all FACTS treatments by Activity Unit - i.e. project
    intersecting_facts_grouped <- intersecting_facts %>%
      group_by(ACTIVITY_UNIT_CN) %>%
      summarise(geometry = st_union(geometry), 
                TREAT_CAT = case_when(length(unique(TREATMENT_CAT)) == 1 ~ first(TREATMENT_CAT),
                                      length(unique(TREATMENT_CAT)) > 1 ~ "Rx & Mech"), 
                YEAR_COMPLETE = max(YEAR, na.rm = T), 
                TREAT_TYPES = paste(unique(TREATMENT_TYPE), collapse = ", "), 
                PLANNED_ACRES = sum(NBR_UNITS_PLANNED, na.rm = TRUE),
                ACCOMPLISHED_ACRES = sum(NBR_UNITS_ACCOMPLISHED, na.rm = TRUE),
                COMPLETE = ifelse(ACCOMPLISHED_ACRES / PLANNED_ACRES < 0.25, 0, 1),
                TREAT_COST = sum(ADJUSTED_TOTAL_COST, na.rm = T))
    
    act_treat_complete <- dplyr::select(as.data.frame(intersecting_facts_grouped), ACTIVITY_UNIT_CN, COMPLETE, TREAT_COST) %>% dplyr::rename(TREAT_ID = ACTIVITY_UNIT_CN)
    
    # Calculate the total acres treated for each treatment project
    intersecting_facts_grouped$Acres_Treated_Treatment <- as.numeric(st_area(intersecting_facts_grouped)/4046.86)
    
    #### Intersect plots with treatments
  
    plots_Treated_Int <- st_intersection(plots, intersecting_facts_grouped)
    
    # calculate acres of treatment in a given plot
    plots_Treated_Int$Acres_Treated_Int <- as.numeric(st_area(plots_Treated_Int)/4046.86)
    
    #### aggregate treatments within a plot to calculate treatment status ("Treated"), treatment type category ("TREAT_CAT"), 
    ####                                        treatment project size ("TREAT_SIZE"), treatment time completed ("YEAR_COMPLETE")
    
    plots_Treated <- plots_Treated_Int %>%
      as.data.frame() %>%
      group_by(Plot_ID) %>%
      summarise(Acres_Treated = sum(Acres_Treated_Int), 
                Grid_Acres = dplyr::first(Grid_Acres), 
                TREAT_CAT = case_when(length(unique(TREAT_CAT)) == 1 ~ first(TREAT_CAT),
                                      length(unique(TREAT_CAT)) > 1 ~ "Rx & Mech"), 
                TREAT_SIZE = sum(Acres_Treated_Treatment), # treatment size is based on the overall size of the treatment project 
                YEAR_COMPLETE = max(YEAR_COMPLETE, na.rm = T), 
                TREAT_ID  = dplyr::first(ACTIVITY_UNIT_CN), 
                TREAT_TYPES = paste(unique(TREAT_TYPES), collapse = ", "))  %>% # 
      mutate(Pct_Treated = Acres_Treated/Grid_Acres,
             Treated = ifelse(Pct_Treated > Tau, 1, 0)) %>% # treatment status is determined by the % threshold Tau 
      dplyr::select(Plot_ID, Pct_Treated, Treated, TREAT_CAT, TREAT_SIZE, YEAR_COMPLETE, TREAT_ID, TREAT_TYPES)
    
    plots_Treated <- merge(plots_Treated, act_treat_complete, by = "TREAT_ID", all.x = T)
    
    plots <- merge(plots, plots_Treated, by = "Plot_ID", all.x = T)
    
    # For plots with no treatment set it equal to 0
    plots[is.na(plots$Treated), "Treated"] <- 0
    plots[is.na(plots$Pct_Treated), "Pct_Treated"] <- 0
    
    #### Create lagged treatment variable - i.e. allow for impact of a the potential spread from a treated plot to a untreated plot
    
    plots_lag_treated <- as.data.frame(plots) %>%
      group_by(direction) %>%   # Group by direction (fire identifier)
      arrange(distance_bin_end) %>%         # Ensure the data is sorted by distance bin
      mutate(Treated_LAG = lag(Treated, n = 1, default = 0)) %>%  # Create the lagged burned variable
      ungroup() %>%
      dplyr::select(direction, distance_bin_end, Treated_LAG) %>%
      arrange(direction)
    
    plots <- merge(plots, plots_lag_treated, 
                   by = c("direction", "distance_bin_end"), all.x = T)
    
    plots <- mutate(plots, distance_bin_end = distance_bin_end/(K*1000)) # normalize the distance bin to kilometers
    
    #### Create indicators for whether a direction is a treated direction or not. Get distance bin of closest treatment in a direction - L_TAU
    
    Direction_L_Taus <- as.data.frame(plots) %>%
      group_by(direction) %>%
      filter(Treated == 1) %>%
      summarise(L_TAU = min(as.numeric(distance_bin_end)),
                Treated_Dir = 1)
    
    plots <- merge(plots, Direction_L_Taus, by = "direction", all.x = T)
    plots <- mutate(plots, AFT_TREAT = ifelse(distance_bin_end > L_TAU & Treated == 0, 1, 0))
    
    plots[is.na(plots$Treated_Dir), "Treated_Dir"] <- 0
    
    ## Indicators of whether a plot is pre or post treatment
    
    plots <- mutate(plots, PRE = ifelse(distance_bin_end < L_TAU, 1, 0),
                    POST = ifelse(distance_bin_end >= L_TAU, 1, 0))
    plots[is.na(plots$PRE), "PRE"] <- 0
    plots[is.na(plots$POST), "POST"] <- 0
    
    plots <- mutate(plots, time_to_treat = ifelse(Treated_Dir == 1, as.numeric(distance_bin_end) - L_TAU, 0))
    
    ## Define Treatment category, size, and time from treatment based on the value for the first treatment to intersect a direction
    
    Treat_Directions <- as.data.frame(plots) %>%
      filter(Treated_Dir == 1 & time_to_treat == 0) %>%
      group_by(direction) %>%
      summarise(TREAT_CAT = dplyr::first(TREAT_CAT), 
                TREAT_SIZE = dplyr::first(TREAT_SIZE), 
                YEAR_COMPLETE = dplyr::first(YEAR_COMPLETE), 
                TREAT_TYPES = dplyr::first(TREAT_TYPES))
    
    plots <- subset(plots, select = -c(TREAT_CAT, TREAT_SIZE, YEAR_COMPLETE, TREAT_TYPES))
    
    plots <- merge(plots, Treat_Directions, by = "direction", all.x = T)
    
    ############   Calculate wildfire outcomes
    
    #### Get Day of Burn & Burn Severity
    
    plots_Cent <- st_centroid(plots)
    
    plots_Cent_reproj <- st_transform(plots_Cent, crs = st_crs(fire))
    
    plots$DAY_BURN <- raster::extract(fire, plots_Cent_reproj) # Extract the day a grid cell burns
    
    #### Intersect wildfire perimeter with plots to get percent burned in a plot
    
    plots_Fire_Int <- st_intersection(plots, mtbs_fire)
    plots_Fire_Int$Fire_Grid_Acres <- as.numeric(st_area(plots_Fire_Int)/4046.86)
    
    plots_Fire_Ints <- plots_Fire_Int %>%
      as.data.frame() %>%
      group_by(Plot_ID) %>%
      summarise(Fire_Grid_Acres = sum(Fire_Grid_Acres), 
                Grid_Acres = dplyr::first(Grid_Acres)) %>%
      mutate(PCT_BURN = Fire_Grid_Acres/Grid_Acres) %>%
      dplyr::select(Plot_ID, PCT_BURN)
    
    plots <- merge(plots, plots_Fire_Ints, by = "Plot_ID", all.x = T)
    plots[is.na(plots$PCT_BURN), "PCT_BURN"] <- 0
    
    ####   Code a grid as burned if the centroid intersects with the fire
    
    plots_Cent_Fire_Int <- st_intersection(plots_Cent, mtbs_fire)
    plots_Cent_Fire_Int <- as.data.frame(plots_Cent_Fire_Int) %>%
      group_by(Plot_ID) %>%
      summarise(BURN = 1)
    
    plots <- merge(plots, plots_Cent_Fire_Int, by = "Plot_ID", all.x = T)  
    plots[is.na(plots$BURN), "BURN"] <- 0
    
    ####   Create lags of whether a plots previous plot in the same direction burns or not ("BURN_LAG").
    
    Fire_Direction_Distance_Burned <- as.data.frame(plots) %>%
      group_by(direction) %>%   # Group by direction (fire identifier)
      arrange(distance_bin_end) %>%         # Ensure the data is sorted by distance bin
      mutate(BURN_LAG = lag(BURN, n = 1, default = 1)) %>%  # Create the lagged burned variable
      ungroup() %>%
      dplyr::select(direction, distance_bin_end, BURN_LAG) %>%
      arrange(direction)
    
    plots <- merge(plots, Fire_Direction_Distance_Burned, 
                   by = c("direction", "distance_bin_end"), all.x = T)
    
    plots$YEAR <- fire_year
    plots$FIRE_SIZE <- fire_size
    
    # Calcualte Year Since Treatment
    
    plots <- mutate(plots, YEARS_SINCE_TREAT = YEAR - YEAR_COMPLETE)
    
    #### For all the plots that don't burn find the raster that it is closest to - allows us to still calculate weather controls on the day of potential burn for unburned plots
    
    plots_No_Burn <- st_transform(plots_Cent, crs(fire))
    
    grid_points <- as(plots_No_Burn, "Spatial")
    
    # Convert the raster to points
    raster_points <- rasterToPoints(fire, spatial = TRUE)
    raster_points <- st_as_sf(raster_points)
    
    # Identify the column name for raster values dynamically
    value_column <- names(raster_points)[sapply(raster_points, is.numeric)]
    
    # Ensure there's only one numeric column that holds raster values
    if(length(value_column) != 1) {
      stop("Unable to uniquely identify the raster value column.")
    }
    
    # Extract the raster values dynamically using the identified column name
    raster_points <- raster_points %>%
      mutate(raster_value = .data[[value_column]])
    
    # Use st_nn to find the nearest neighbor raster point for each grid centroid
    nearest_indices <- st_nn(st_as_sf(plots_No_Burn), raster_points, returnDist = TRUE)
    
    # Extract the nearest raster values
    nearest_values <- sapply(nearest_indices$nn, function(idx) raster_points$raster_value[idx])
    
    # Add the raster values to the original plots_Cent data frame
    plots_No_Burn <- plots_No_Burn %>%
      mutate(DAY_BURN_CLOSEST = as.numeric(nearest_values))
    
    plots <- merge(plots, dplyr::select(as.data.frame(plots_No_Burn), Plot_ID, DAY_BURN_CLOSEST), by = "Plot_ID", all.x = T)
    
    # Define date burned as the date where the closest grid burns
    plots <- mutate(plots, DATE_BURNED = as.Date(paste0(YEAR, "-", DAY_BURN_CLOSEST), format = "%Y-%j")) 
    
    plots_reproj <- st_transform(plots, crs = st_crs(mtbsBS_year))
    
    
    ###### Extract mean ordinal MTBS Burn Severity ("BURN_SEV") and the percent of plot burned at moderate-high and high burn severity 
    ######                                                            ("BURN_SEV_MED_HIGH_PCT" and "BURN_SEV_HIGH_PCT") for each grid.
    
    plots$BURN_SEV <- as.numeric(exact_extract(mtbsBS_year, plots_reproj, "mean", progress = F))
    
    plots <- mutate(plots, BURN_SEV0 = ifelse(is.na(BURN_SEV) == T, 0, BURN_SEV)) # BURN_SEV0 codes non-burned plots as zeros
    
    #### Calculate Moderate-High Burn Severity
    
    burn_sev_34_pct <- function(df) {
      valid <- !is.na(df$value)
      num <- sum(df$coverage_fraction[valid & df$value %in% c(3, 4)])
      denom <- sum(df$coverage_fraction[valid])
      if (denom == 0) return(NA) else return(num / denom)
    }
    
    
    plots$BURN_SEV_MED_HIGH_PCT <- exact_extract(
      mtbsBS_year,
      plots_reproj,
      fun = burn_sev_34_pct,
      summarize_df = TRUE,
      progress = FALSE
    )
    
    plots[is.na(plots$BURN_SEV_MED_HIGH_PCT), "BURN_SEV_MED_HIGH_PCT"] <- 0
    
    #### Calculate High Burn Severity
    
    burn_sev_4_pct <- function(df) {
      valid <- !is.na(df$value)
      num <- sum(df$coverage_fraction[valid & df$value %in% c(4)])
      denom <- sum(df$coverage_fraction[valid])
      if (denom == 0) return(NA) else return(num / denom)
    }
    
    
    plots$BURN_SEV_HIGH_PCT <- exact_extract(
      mtbsBS_year,
      plots_reproj,
      fun = burn_sev_4_pct,
      summarize_df = TRUE,
      progress = FALSE
    )
    
    plots[is.na(plots$BURN_SEV_HIGH_PCT), "BURN_SEV_HIGH_PCT"] <- 0
    
    
    plots$FIRE_ID <- fire_id
    plots$FIRE_NAME <- fire_name
    
    plots <- plots[, c(setdiff(names(plots), "geometry"), "geometry")] # make sure geometry is last column
    plots <- dplyr::rename(plots, distance_bin = distance_bin_end)
    
    #### Create an indicator if the fire in a treated direction ever hits the treatment ("HIT_TREAT").
    
    Treat_Dir_Hit <- as.data.frame(plots) %>%
      filter(Treated_Dir == 1) %>%
      group_by(direction) %>%
      summarise(L_TAU = dplyr::first(L_TAU), 
                Max_Dist_Burn = max(distance_bin[BURN == 1])) %>%
      mutate(HIT_TREAT = ifelse(L_TAU <= Max_Dist_Burn + 1, 1, 0)) %>%
      dplyr::select(direction, HIT_TREAT)
    
    plots <- merge(plots, Treat_Dir_Hit, by = "direction", all.x = T)
    
    plots[is.na(plots$HIT_TREAT), "HIT_TREAT"] <- 0
    
    #### For each treatment-direction get the distance bin it first hits the treatment - "DIST_TREAT"
    
    Treat_Dist_Hit <- as.data.frame(plots) %>%
      filter(Treated_Dir == 1) %>%
      group_by(direction, TREAT_ID) %>%
      summarise(DIST_TREAT = min(distance_bin[Treated == 1])) %>%
      dplyr::select(direction,TREAT_ID, DIST_TREAT)
    
    Treat_Dist_Hit <- Treat_Dist_Hit[!is.na(Treat_Dist_Hit$TREAT_ID), ]
    
    plots <- merge(plots, Treat_Dist_Hit, by = c("direction", "TREAT_ID"), all.x = T)
    
    plots[is.na(plots$DIST_TREAT), "DIST_TREAT"] <- 0
    
    ## Get the distance bins for each treatment that intersects with a fire (conditional not yet being extinguished - BURN_LAG == 1) in a given direction
    
    Direction_Distance_Treatments <- as.data.frame(plots) %>%
      filter(!is.na(TREAT_ID)) %>%
      group_by(direction, TREAT_ID) %>%
      filter(Treated == 1) %>%
      summarise(Min_Dist = min(as.numeric(distance_bin)),
                BURN_LAG = max(BURN_LAG), 
                HIT_TREAT = max(HIT_TREAT)) %>%
      filter(BURN_LAG != 0) %>%
      group_by(direction) %>%
      summarise(DIST_HIT_TREATS = paste(unique(Min_Dist), collapse = ", "))
    
    plots <- merge(plots, Direction_Distance_Treatments, by = "direction", all.x = T)
    
    #### For each Grid get an indicator if a grid cell has previously burned in the last 10 years
    
    mtbs_previous_10yr <- filter(mtbs, fire_year > YEAR_MTBS) %>% filter(fire_year - YEAR_MTBS <= 10)
    plots_Burned_Int <- st_intersection(plots, mtbs_previous_10yr)
    plots_Burned_Int$Acres_Burned_Int <- as.numeric(st_area(plots_Burned_Int)/4046.86)
    
    plots_Burned <- plots_Burned_Int %>%
      as.data.frame() %>%
      group_by(Plot_ID) %>%
      summarise(Previous_Acres_Burned = sum(Acres_Burned_Int), 
                Grid_Acres = dplyr::first(Grid_Acres))  %>%
      mutate(Pct_Burned_10Y = Previous_Acres_Burned/Grid_Acres,
             PREV_BURN_10Y = ifelse(Pct_Burned_10Y > 0, 1, 0)) %>%
      dplyr::select(Plot_ID, Pct_Burned_10Y, PREV_BURN_10Y)
    
    plots <- merge(plots, plots_Burned, by = "Plot_ID", all.x = T)
    
    plots[is.na(plots$Pct_Burned_10Y), "Pct_Burned_10Y"] <- 0
    plots[is.na(plots$PREV_BURN_10Y), "PREV_BURN_10Y"] <- 0
    
    #### Create fire-direction, fire-distance, fire-direction-distance FEs, & Burn Right Away Indicator (i.e. if fire intersects with treatment right away)
    
    plots <- mutate(plots, direction_distance_fire_FE = as.factor(paste(direction, distance_bin, FIRE_ID, sep = "-")),
                    direction_fire_FE = paste(direction, FIRE_ID, sep = "-"),
                    distance_fire_FE = paste(distance_bin, FIRE_ID, sep = "-"), 
                    Burn_Right_Away = ifelse(is.na(L_TAU) == F & L_TAU == 1, 1, 0))
    
    #### Extract MTT outputs
    
    plots_reproj <- st_transform(plots, crs = st_crs(MTT))
  
    plots$MTT <- as.numeric(exact_extract(MTT, plots_reproj, "mean", progress = F))
    
    Fire_Direction_Distance_MTT <- as.data.frame(plots) %>%
      group_by(direction) %>%   # Group by direction (fire identifier)
      arrange(distance_bin) %>%         # Ensure the data is sorted by distance bin
      mutate(MTT_LAG = lag(MTT, n = 1, default = NA)) %>%  # Create the lagged burned variable
      ungroup() %>%
      dplyr::select(direction, distance_bin, MTT_LAG) %>%
      arrange(direction, distance_bin)
    
    plots <- merge(plots, Fire_Direction_Distance_MTT, 
                   by = c("direction", "distance_bin"), all.x = T)  
    
    plots <- mutate(plots, DELTA_MTT = ifelse(distance_bin == 1, MTT, MTT - MTT_LAG))
    
    ## Delta MTT Missing Indicator - missing if there is no MTT in previous cell do to missing fuel or if arrival time is lower than in focal cell
    
    plots <- mutate(plots, NA_DELTA_MTT = ifelse(is.na(DELTA_MTT) == T | DELTA_MTT < 0, 1, 0))
    
    #### Extract Fireline Intensity from MTT
    
    plots_reproj <- st_transform(plots, crs = st_crs(MTT_INT))
    
    plots$FIRE_INTENSITY <- as.numeric(exact_extract(MTT_INT, plots_reproj, "mean", progress = F))
    plots <- mutate(plots, LOG_FIRE_INTENSITY = log(FIRE_INTENSITY))
    
    
    #### Calculate Distance to Suppression lines
    
    QAQC_fire <- filter(QAQC, MTBS_id == fire_id)
    
    if (nrow(QAQC_fire) == 0){
      
      plots <- mutate(plots, SUP_LINES_PRESENT = 0)
      plots$SUP_LINE_DIST <- NA
      plots$SUP_LINE <- NA
      plots$LOG_SUP_LINE_DIST <- NA
      plots$SUP_LINE_LEN <- NA
      plots$SUP_LINE_INT <- NA
    }
    
    if (nrow(QAQC_fire) != 0){
      
      plots <- mutate(plots, SUP_LINES_PRESENT = 1)
      plots$SUP_LINE_DIST <- calculate_distance(plots, QAQC_fire)
      plots <- mutate(plots, SUP_LINE = ifelse(SUP_LINE_DIST == 0, 1, 0),
                      LOG_SUP_LINE_DIST = log(1 + SUP_LINE_DIST))
      
      
      ## Calculate the length of suppression lines in each grid
      
      # Find the intersections of lines with each grid cell
      lines_in_plots <- st_intersection(QAQC_fire, plots)
      
      # Calculate the length of each intersected line segment
      lines_in_plots <- lines_in_plots %>%
        mutate(segment_length = st_length(geometry))
      
      # Sum the lengths for each grid cell
      grid_line_lengths <- as.data.frame(lines_in_plots) %>%
        group_by(Plot_ID) %>%
        summarise(SUP_LINE_LEN = as.numeric(sum(segment_length)))
      
      plots <- merge(plots, grid_line_lengths, by = "Plot_ID", all.x = T)
      
      plots[is.na(plots$SUP_LINE_LEN), "SUP_LINE_LEN"] <- 0
      
      plots <- mutate(plots, SUP_LINE_INT = SUP_LINE_LEN/Grid_Acres)
      
    }
    
    #### Create unique identifier (ID) defined as the combination of the plot and fire ID 
    
    plots <- mutate(plots, ID = paste(Plot_ID, FIRE_ID, sep = "-"))
    plots <- subset(plots, select = -c(distance_bin_start, Plot_ID, YEAR_COMPLETE))
    
    plots <- plots %>%
      dplyr::select(ID, everything())
    
    colnames(plots)
    colnames(panel)
    
    panel <- rbind(panel, plots)
  }
  
  toc()
  
  #### Merge smoke data to panel
  
  panel <- merge(panel, smoke_data, by = "FIRE_ID", all.x = T)
  
  ## Create plot level CO2 and PM2.5 per acre emissions estimates by assuming emissions are proportional to fire size 
  
  panel <- mutate(panel, FIRE_PM25_PER_ACRE = FIRE_PM25/FIRE_SIZE, 
                             FIRE_CO2_PER_ACRE = FIRE_CO2/FIRE_SIZE)
  
  panel <- merge(panel, smoke_costs, by = "FIRE_ID", all.x = T)
  
  ## Similarly for per acre death and earning estimates
  
  panel <- mutate(panel, Deaths_PER_ACRE = Deaths/FIRE_SIZE,
                             DeathCost_PER_ACRE = DeathCost/FIRE_SIZE, 
                             EarningsLoss_PER_ACRE = EarningsLoss/FIRE_SIZE)
  
  ################ Step 4. Extract other plot characteristics from the full panel that need parallel processing
  
  pacman::p_load(dplyr,sf, tmap, magrittr, rnaturalearth, rnaturalearthdata, ggplot2, maps, lwgeom, raster, stars, haven, stargazer, quantmod, lubridate, tidyr, ggpubr,
                exactextractr, tictoc, terra, gtools, here, fixest, modelsummary, readr, rdrobust, prism, parallel,tmaptools,
                 OpenStreetMap, maptiles, gifski, stringr, tidycensus)
  
  
  OLS_Grids <- panel
  
  ######## Calculate Observable Characteristics of Grids
  
  #### The function "extract_gridMET_data" extracts the wind speed, wind direction, ERC, and 1000 hr fuel moisture based on the day of burn for a plot and its centroid
  
  extract_gridMET_data <- function(OLS_Grids){
    
    df <- data.frame(matrix(ncol = 6, nrow = 0))
    
    colnames(df) <- c('ID', 'DATE_BURNED', 'WindSpeed', 'WindDirection', 'ERC', 'FM1000')
    
    OLS_Grids_Cent <- st_centroid(OLS_Grids)
    
    for (yr in unique(OLS_Grids$YEAR)){
      
      df_temp <- data.frame(matrix(ncol = 3, nrow = 0))
      colnames(df_temp) <- c('ID', 'Row_ID', 'DATE_BURNED')
      
      OLS_Grids_yr <- filter(OLS_Grids_Cent, YEAR == yr)
      
      grids_ID <- OLS_Grids_yr$ID
      
      df_temp[1:length(grids_ID),1] <- grids_ID
      df_temp$Row_ID <- 1:nrow(OLS_Grids_yr)
      
      burn_dates <- as.Date(OLS_Grids_yr$DATE_BURNED, "%m/%d/%Y")
      df_temp$DATE_BURNED <- burn_dates
      
      vs_yr <- rast(here("data", "raw", "gridMET", paste("vs_", yr, ".nc", sep = "")))
      erc_yr <- rast(here("data", "raw", "gridMET", paste("erc_", yr, ".nc", sep = "")))
      th_yr <- rast(here("data", "raw", "gridMET", paste("th_", yr, ".nc", sep = "")))
      fm1000_yr <- rast(here("data", "raw", "gridMET", paste("fm1000_", yr, ".nc", sep = "")))
      
      OLS_Grids_yr_reproj <- st_transform(OLS_Grids_yr, crs = st_crs(vs_yr))
      
      Grids_WindSpeed <- raster::extract(vs_yr, OLS_Grids_yr_reproj) %>%
        pivot_longer(cols = -c(ID), names_to = "date_type") %>%
        dplyr::rename(Row_ID = ID) %>%
        mutate(num = as.numeric(str_sub(date_type, -5, -1))) %>%
        mutate(date = num + lubridate::ymd("1900-01-01")) %>%
        filter(date %in% burn_dates) %>%
        arrange(Row_ID) %>%
        dplyr::rename(WindSpeed = value, DATE_BURNED = date) %>%
        dplyr::select(Row_ID, DATE_BURNED, WindSpeed)
      
      Grids_WindDirection <- raster::extract(th_yr, OLS_Grids_yr_reproj) %>%
        pivot_longer(cols = -c(ID), names_to = "date_type") %>%
        dplyr::rename(Row_ID = ID) %>%
        mutate(num = as.numeric(str_sub(date_type, -5, -1))) %>%
        mutate(date = num + lubridate::ymd("1900-01-01")) %>%
        filter(date %in% burn_dates) %>%
        arrange(Row_ID) %>%
        dplyr::rename(WindDirection = value, DATE_BURNED = date) %>%
        dplyr::select(Row_ID, DATE_BURNED, WindDirection)
      
      df_temp <- merge(df_temp, Grids_WindSpeed, by = c("Row_ID", "DATE_BURNED"), all.x = T)
      df_temp <- merge(df_temp, Grids_WindDirection, by = c("Row_ID", "DATE_BURNED"), all.x = T)
      
      
      Grids_ERC <- raster::extract(erc_yr, OLS_Grids_yr_reproj) %>%
        pivot_longer(cols = -c(ID), names_to = "date_type") %>%
        dplyr::rename(Row_ID = ID) %>%
        mutate(num = as.numeric(str_sub(date_type, -5, -1))) %>%
        mutate(date = num + lubridate::ymd("1900-01-01")) %>%
        filter(date %in% burn_dates) %>%
        arrange(Row_ID) %>%
        dplyr::rename(ERC = value, DATE_BURNED = date) %>%
        dplyr::select(Row_ID, DATE_BURNED, ERC)
      
      df_temp <- merge(df_temp, Grids_ERC, by = c("Row_ID", "DATE_BURNED"), all.x = T)
      
      Grids_fm1000 <- raster::extract(fm1000_yr, OLS_Grids_yr_reproj) %>%
        pivot_longer(cols = -c(ID), names_to = "date_type") %>%
        dplyr::rename(Row_ID = ID) %>%
        mutate(num = as.numeric(str_sub(date_type, -5, -1))) %>%
        mutate(date = num + lubridate::ymd("1900-01-01")) %>%
        filter(date %in% burn_dates) %>%
        arrange(Row_ID) %>%
        dplyr::rename(FM1000 = value, DATE_BURNED = date) %>%
        dplyr::select(Row_ID, DATE_BURNED, FM1000)
      
      df_temp <- merge(df_temp, Grids_fm1000, by = c("Row_ID", "DATE_BURNED"), all.x = T)
      
      df_temp <- dplyr::select(df_temp, ID, DATE_BURNED, WindSpeed, WindDirection, ERC, FM1000)
      df <- rbind(df, df_temp)
    }
    
    return(df)
  }
  
  
  ## Set up number of cores for parallel processing
  
  num_cores <- 7
  
  num_in_group <- floor(nrow(OLS_Grids) / num_cores)
  
  #--- assign group id to polygons ---#
  OLS_Grids <- OLS_Grids %>%
    mutate(
      #--- create grid id ---#
      grid_id = 1:nrow(.),
      #--- assign group id  ---#
      group_id = grid_id %/% num_in_group + 1
    )
  
  #### i) Weather Variables
  
  tic()
  
  #--- parallelized processing by group ---#
  temp <- mclapply(
    1:(num_cores + 1),
    function(i) extract_gridMET_data(filter(OLS_Grids, group_id == i)),
    mc.cores = (num_cores + 1)
  )
  toc()
  
  OLS_Grids <- merge(OLS_Grids, dplyr::select(bind_rows(temp), ID, WindSpeed, WindDirection, ERC, FM1000), by = "ID", all.x = T)
  
  
  #### Wind Difference - This is the Cosine of the difference in the bearing of a grids direction minus the dominant wind direction on the day the previous cells burned
  
  Fire_Direction_Wind_Direction <- as.data.frame(OLS_Grids) %>%
    group_by(direction_fire_FE) %>%   # Group by direction (fire identifier)
    arrange(distance_bin) %>%         # Ensure the data is sorted by distance bin
    mutate(Wind_Dir_LAG = lag(WindDirection, n = 1, default = NA)) %>%  # Create the lagged burned variable
    ungroup() %>%
    dplyr::select(direction_fire_FE, distance_bin, Wind_Dir_LAG) %>%
    arrange(direction_fire_FE, distance_bin)
  
  OLS_Grids <- merge(OLS_Grids, Fire_Direction_Wind_Direction, 
                     by = c("direction_fire_FE", "distance_bin"), all.x = T) 
  
  
  OLS_Grids <- mutate(OLS_Grids, WIND_DIFF = ifelse(distance_bin != 1, cos((bearing - Wind_Dir_LAG)*pi/180), 
                                                    cos((bearing - WindDirection)*pi/180))
  )
  
  #### ii) Extract topographic variables from LANDFIRE
  
  tic()
  
  OLS_Grids_reproj <- st_transform(OLS_Grids, crs = st_crs(slope))
  
  #--- parallelized processing by group ---#
  temp <- mclapply(
    1:(num_cores + 1),
    function(i) exact_extract(top_landfire, filter(OLS_Grids_reproj, group_id == i), "mean", progress = F) %>%
      as.data.frame() %>%
      dplyr::rename(Slope = mean.LC20_SlpD_220,
                    Aspect = mean.LC20_Asp_220,
                    Elev = mean.LC20_Elev_220,
                    TRI = mean.tri),
    mc.cores = (num_cores + 1)
  )
  OLS_Grids <- cbind(OLS_Grids, bind_rows(temp))
  
  rm(slope, elev, aspect, TRI) # remove loaded topography datasets
  
  OLS_Grids <- mutate(OLS_Grids, Aspect_Class = case_when(is.na(Aspect) ~ "No Aspect",
                                                          Aspect == -1 ~ "Flat",
                                                          Aspect >= 337.5 | Aspect < 22.5 ~ "North",
                                                          Aspect >= 22.5 & Aspect < 67.5 ~ "Northeast",
                                                          Aspect >= 67.5 & Aspect < 112.5 ~ "East",
                                                          Aspect >= 112.5 & Aspect < 157.5 ~ "Southeast",
                                                          Aspect >= 157.5 & Aspect < 202.5 ~ "South",
                                                          Aspect >= 202.5 & Aspect < 247.5 ~ "Southwest",
                                                          Aspect >= 247.5 & Aspect < 292.5 ~ "West",
                                                          Aspect >= 292.5 & Aspect < 337.5 ~ "Northwest",
                                                          TRUE ~ "Unknown")) %>%
    mutate(SOUTH_FACING = ifelse(Aspect_Class == "Southeast" | Aspect_Class == "South" | Aspect_Class == "Southwest", 1, 0))
  
  toc()
  
  #### iii) Extract fire regime characteristics (MFRI) and fuel model type from LANDFIRE 
  
  # MFRI
  
  OLS_Grids_reproj <- st_transform(OLS_Grids, crs = st_crs(MFRI))
  
  temp <- mclapply(
    1:(num_cores + 1),
    function(i) as.numeric(exact_extract(MFRI, filter(OLS_Grids_reproj, group_id == i), "mean", progress = F)),
    mc.cores = (num_cores + 1)
  )
  OLS_Grids$MFRI <- unlist(temp)
  
  rm(MFRI)
  
  # Fuel Model Type
  
  OLS_Grids_Centroids <- st_centroid(subset(OLS_Grids, select = c(ID,group_id, geometry)))
  
  OLS_Grids_Centroids_reproj <- st_transform(OLS_Grids_Centroids, crs = st_crs(FBFM13))
  
  tic()
  
  temp <- mclapply(
    1:(num_cores + 1),
    function(i) as.numeric(raster::extract(FBFM13, filter(OLS_Grids_Centroids_reproj, group_id == i))),
    mc.cores = (num_cores + 1)
  )
  OLS_Grids$FBFM13 <- unlist(temp)
  
  rm(FBFM13)
  
  # Will categorize each grid into either "Grass", "Shrub", "Timber", "Slash", or "Other" Fuel Types based on - https://www.fs.usda.gov/rm/pubs/rmrs_gtr175/rmrs_gtr175_367_396.pdf
  
  OLS_Grids <- mutate(OLS_Grids,
                      FuelType_2001 = case_when(FBFM13 == 1 | FBFM13 == 2 | FBFM13 == 3 ~ "Grass",
                                                FBFM13 == 4 | FBFM13 == 5 | FBFM13 == 6 | FBFM13 == 7 ~ "Shrub",
                                                FBFM13 == 8 | FBFM13 == 9 | FBFM13 == 10 ~ "Timber",
                                                FBFM13 == 11 | FBFM13 == 12 | FBFM13 == 13 ~  "Slash",
                                                FBFM13 == 91 | FBFM13 == 92 | FBFM13 == 93 | FBFM13 == 98 | FBFM13 == 99 | FBFM13 == -9999 ~ "Other"))
  
  
  toc()
  
  
  #### iv) Calculate Distance to USFS Road, US Highway, and WUI
  
  ## Distance to USFS Road
  tic()
  
  temp <- mclapply(
    1:(num_cores + 1),
    function(i) calculate_distance(filter(OLS_Grids, group_id == i), USFS_Roads),
    mc.cores = (num_cores + 1)
  )
  OLS_Grids$Distance_FS_Road <- unlist(temp)
  
  OLS_Grids <- mutate(OLS_Grids, Contains_USFS_Road = ifelse(Distance_FS_Road == 0, 1, 0))
  toc()
  
  ## Distance to US Highway
  
  tic()
  
  temp <- mclapply(
    1:(num_cores + 1),
    function(i) calculate_distance(filter(OLS_Grids, group_id == i), US_Highways),
    mc.cores = (num_cores + 1)
  )
  OLS_Grids$Distance_US_Highway <- unlist(temp)
  
  toc()
  
  OLS_Grids <- mutate(OLS_Grids, ROAD = ifelse(Distance_US_Highway == 0 | Distance_FS_Road == 0, 1, 0))
  
  
  ## Distance to WUI
  
  tic()
  
  temp <- mclapply(
    1:(num_cores + 1),
    function(i) calculate_distance(filter(OLS_Grids, group_id == i), Western_WUI),
    mc.cores = (num_cores + 1)
  )
  OLS_Grids$Distance_WUI <- unlist(temp)
  
  toc()
  
  OLS_Grids <- mutate(OLS_Grids, WUI = ifelse(Distance_WUI == 0, 1, 0))
  
  
  #### v) Calculate No. Strucutures & Housing Units inside of a grid
  
  # Building Count
  
  tic()
  
  
  temp <- mclapply(
    1:(num_cores + 1),
    function(i) as.numeric(exact_extract(Buildings, filter(OLS_Grids_reproj, group_id == i), 'sum', progress = F)),
    mc.cores = (num_cores + 1)
  )
  toc()
  
  OLS_Grids$Struc_Count <- unlist(temp)
  
  rm(Buildings)
  
  
  # Housing Count
  
  OLS_Grids_reproj <- st_transform(OLS_Grids, crs = st_crs(HUCount))
  
  tic()
  
  temp <- mclapply(
    1:(num_cores + 1),
    function(i) as.numeric(exact_extract(HUCount, filter(OLS_Grids_reproj, group_id == i), 'sum', progress = F)),
    mc.cores = (num_cores + 1)
  )
  toc()
  
  OLS_Grids$HU_Count <- unlist(temp)
  
  rm(HUCount)
  
  
  ## For grids with Non-zero buildings or structures calculate the median housing value from ACS
  
  OLS_Grids_struc <- filter(OLS_Grids, HU_Count > 0 | Struc_Count > 0) %>% st_centroid()
  
  ## Calculate housing value corresponding to year of burn. For observations that have NA default to housing value of 2016-2020 survey
  
  Grids_ACS_Int <- st_intersection(OLS_Grids_struc, ACS_sf) %>%
    as.data.frame() %>%
    mutate(across(starts_with("MED_HOUSE_VAL_"), as.numeric)) %>%  # Ensure numeric values
    mutate(MED_HOUSE_VAL = case_when(
      YEAR == 2017 ~ MED_HOUSE_VAL_2017,
      YEAR == 2018 ~ MED_HOUSE_VAL_2018,
      YEAR == 2019 ~ MED_HOUSE_VAL_2019,
      YEAR == 2020 ~ MED_HOUSE_VAL_2020,
      YEAR == 2021 ~ MED_HOUSE_VAL_2021,
      YEAR == 2022 | YEAR == 2023 ~ MED_HOUSE_VAL_2022
    )) %>%
    mutate(MED_HOUSE_VAL = ifelse(is.na(MED_HOUSE_VAL) == T, MED_HOUSE_VAL_2020, MED_HOUSE_VAL)) %>%
    dplyr::select(ID, MED_HOUSE_VAL)
  
  OLS_Grids <- merge(OLS_Grids, Grids_ACS_Int, by = "ID", all.x = T)
  
  OLS_Grids <- mutate(OLS_Grids, TOT_HOUSE_VAL = MED_HOUSE_VAL*HU_Count, TOT_STRUC_VAL = MED_HOUSE_VAL*Struc_Count)
  
  #### vi) Fire Suppression Effort Variables
  
  ## Get Distance to LAT Drops
  
  OLS_Grids_final <- OLS_Grids[0,]
    
  for (yr in seq(2017, 2023, 1)){
    
    print(yr)
    
    OLS_Grids_yr <- filter(OLS_Grids, YEAR == yr)
    OLS_Grids_yr$DIST_LAT <- calculate_distance(OLS_Grids_yr, filter(LAT, YEAR == yr))
    OLS_Grids_final <- rbind(OLS_Grids_final, OLS_Grids_yr)
    
  }
  
  OLS_Grids_final <- mutate(OLS_Grids_final, LAT = ifelse(DIST_LAT == 0 & is.na(DIST_LAT) == F, 1, 0))
  OLS_Grids_final <- mutate(OLS_Grids_final, LOG_DIST_LAT = log(DIST_LAT + 1))
  
  
  
  # Find the intersections of lines with each grid cell
  lines_in_grids <- st_intersection(LAT, filter(OLS_Grids_final, LAT == 1))
  
  # Calculate the length of each intersected line segment
  lines_in_grids <- lines_in_grids %>%
    mutate(segment_length = st_length(geometry))
  
  # Sum the lengths for each grid cell
  grid_line_lengths <- as.data.frame(lines_in_grids) %>%
    group_by(ID) %>%
    summarise(LAT_LINE_LEN = as.numeric(sum(segment_length)))
  
  OLS_Grids_final <- merge(OLS_Grids_final, grid_line_lengths, by = "ID", all.x = T)
  
  OLS_Grids_final[is.na(OLS_Grids_final$LAT_LINE_LEN), "LAT_LINE_LEN"] <- 0
  
  OLS_Grids_final <- mutate(OLS_Grids_final, LAT_LINE_INT = LAT_LINE_LEN/Grid_Acres)
    
  plot_panel_sf <- OLS_Grids_final %>%
    dplyr::select(ID)
  
  ## Save plots as dataframe
  plot_panel_df <- as.data.frame(OLS_Grids_final) %>% subset(select = -c(geometry, grid_id, group_id))
  
  ## Return both shapefile and dataframe
  
  return(list(
    plot_panel_df = plot_panel_df,
    plot_panel_sf = plot_panel_sf
  ))
  
}

#### Run baseline sample with no. directions L = 24, distance bin distance K = 0.5 km, % treatment threshold = 50%, and imputed ignition points from day 1 centroids 

plot_panel_output <- create_SpatialDiD_plot_panel(L = 24, K = 0.5, Tau = 0.5, CENT = TRUE)

plot_panel_df_L24_K05 <- plot_panel_output$plot_panel_df
plot_panel_df_L24_K05_sf <- plot_panel_output$plot_panel_sf

#### Save the dataset

write_csv(plot_panel_df_L24_K05, here("data", "temp", "SpatialDiD_Grids_L24_K05.csv")) # Put the main sample in the "temp" directory to make sure to do smoke imputation in "DescriptiveStats.R" before analysis

toc()

#### Save the shapefile of the baseline plots

dir.create(here("data", "intermediate", "SpatialDiD_Grids_L24_K05"))

st_write(plot_panel_df_L24_K05_sf, here("data", "intermediate", "SpatialDiD_Grids_L24_K05", "SpatialDiD_Grids_L24_K05.shp"), append = FALSE)

######## Optional: create panels to be used in robustness checks. Uncomment below if you want to run these

# #### Robustness Check 1 - Changing No. Directions
# 
# plot_panel_df_L36_K05 <- create_SpatialDiD_plot_panel(L = 36, K = 0.50, Tau = 0.5, CENT = TRUE)$plot_panel_df
# plot_panel_df_L18_K05 <- create_SpatialDiD_plot_panel(L = 18, K = 0.50, Tau = 0.5, CENT = TRUE)$plot_panel_df
# plot_panel_df_L12_K05 <- create_SpatialDiD_plot_panel(L = 12, K = 0.50, Tau = 0.5, CENT = TRUE)$plot_panel_df
# 
# write_csv(plot_panel_df_L36_K05, here("data", "intermediate", "SpatialDiD_Grids_L36_K05.csv"))
# write_csv(plot_panel_df_L18_K05, here("data", "intermediate", "SpatialDiD_Grids_L18_K05_40p.csv")) # treatment thresholds change slightly to account for the change in size of grids when lowering the number of directions.
# write_csv(plot_panel_df_L12_K05, here("data", "intermediate", "SpatialDiD_Grids_L12_K05_30p.csv"))
# 
# #### Robustness Check 2 - Reported ignitions
# 
# plot_panel_df_L24_K05_V2 <- create_SpatialDiD_plot_panel(L = 24, K = 0.50, Tau = 0.5, CENT = FALSE)$plot_panel_df
# 
# write_csv(plot_panel_df_L24_K05_V2, here("data", "intermediate", "SpatialDiD_Grids_L24_K05_V2.csv"))
# 
# #### Robustness Check 3 - Alternative Treatment Thresholds
# 
# plot_panel_df_L24_K05_100p <- create_SpatialDiD_plot_panel(L = 24, K = 0.50, Tau = 0.99, CENT = TRUE)$plot_panel_df
# plot_panel_df_L24_K05_75p <- create_SpatialDiD_plot_panel(L = 24, K = 0.50, Tau = 0.75, CENT = TRUE)$plot_panel_df
# plot_panel_df_L24_K05_25p <- create_SpatialDiD_plot_panel(L = 24, K = 0.50, Tau = 0.25, CENT = TRUE)$plot_panel_df
# plot_panel_df_L24_K05_zero <- create_SpatialDiD_plot_panel(L = 24, K = 0.50, Tau = 0, CENT = TRUE)$plot_panel_df
# 
# write_csv(plot_panel_df_L24_K05_100p, here("data", "intermediate", "SpatialDiD_Grids_L24_K05_100p.csv"))
# write_csv(plot_panel_df_L24_K05_75p, here("data", "intermediate", "SpatialDiD_Grids_L24_K05_75p.csv"))
# write_csv(plot_panel_df_L24_K05_25p, here("data", "intermediate", "SpatialDiD_Grids_L24_K05_25p.csv"))
# write_csv(plot_panel_df_L24_K05_zero, here("data", "intermediate", "SpatialDiD_Grids_L24_K05_zero.csv"))







