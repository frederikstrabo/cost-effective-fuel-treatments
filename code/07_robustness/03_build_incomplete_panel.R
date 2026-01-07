##########################################    03_build_incomplete_panel  ##########################################

################# Purpose: Creates sample of grids to be used for Spatial DiD Analysis

rm(list=ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr,sf, tmap, magrittr, rnaturalearth, rnaturalearthdata, ggplot2, maps, lwgeom, rgeos, raster, stars, haven, stargazer, quantmod, lubridate, tidyr, ggpubr, 
               rgdal, exactextractr, tictoc, terra, gtools, here, fixest, modelsummary, readr, rdrobust, prism, parallel,tmaptools, 
               OpenStreetMap, maptiles, gifski, purrr, nngeo, stringr, geosphere)


# Set Path
here::i_am("code/07_robustness/03_build_incomplete_panel.R")

# Load Functions 

source(here("code", "functions","tidy_facts.R"))
source(here("code", "functions","calculate_distance.R"))
st_erase = function(x, y) st_difference(x, st_make_valid(st_union(st_combine(y)))) # Taken from here https://r-spatial.github.io/sf/reference/geos_binary_ops.html

########    Load in Datasets    ########

#### FACTS - Forest Service Fuel Treatments

facts <- st_read(here("data", "raw", "FACTS", "S_USA.Activity_HazFuelTrt_PL.shp"))
facts <- tidy_facts(facts, inf_yr = 2020, crs = 5070)
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

# Load in USFS NF

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

#### Load in list of fires that intersect with FTs

fires_int_df <- read_csv(here("data", "intermediate", "FACTS_Parks_Fire_List.csv"))

#### Load WFEIS Fire Smoke Emissions

smoke_data <- read_csv(here("data", "intermediate", "WFEIS", "WFEIS_Emissions.csv"))

###### Step 1. Construct set of grids to use for Spatial RD Analysis

# Set up intial sf object to save all grids in

SpatialOLS_Grids <- data.frame(matrix(ncol = 50, nrow = 0), geometry = st_sfc())

colnames(SpatialOLS_Grids) <- c("ID", "direction", "distance_bin", "TREAT_ID", "bearing", "USFS_NF", "Wilderness", "Grid_Acres", "Pct_Treated", "Treated", "COMPLETE",
                                "Treated_LAG", "L_TAU", "Treated_Dir", "AFT_TREAT", "PRE", "POST", "time_to_treat", "TREAT_CAT", "TREAT_SIZE", "TREAT_TYPES", "DAY_BURN", 
                                "PCT_BURN", "BURN", "BURN_LAG", "YEAR", "FIRE_SIZE", "YEARS_SINCE_TREAT", "DAY_BURN_CLOSEST", "DATE_BURNED", "BURN_SEV", "BURN_SEV0", 
                                "BURN_SEV_MED_HIGH_PCT", 
                                "FIRE_ID", "FIRE_NAME", "HIT_TREAT","DIST_TREAT", "DIST_HIT_TREATS", "Pct_Burned_10Y", "PREV_BURN_10Y", 
                                "direction_distance_fire_FE", "direction_fire_FE", "distance_fire_FE", "Burn_Right_Away", "MTT", "MTT_LAG", "DELTA_MTT", "NA_DELTA_MTT",
                                "FIRE_INTENSITY", "LOG_FIRE_INTENSITY", "geometry")


SpatialOLS_Grids <- st_as_sf(SpatialOLS_Grids, crs = 5070)

for (i in 1:nrow(fires_int_df)){
  
  print(i)
  
  n = i
  
  fire_id <- as.character(fires_int_df[i, 1])
  ras_id <- as.numeric(fires_int_df[i, 2])
  
  fire <- raster_list[[ras_id]]
  mtbs_fire <- filter(mtbs, FIRE_ID == fire_id)
  fire_name <- mtbs_fire$FIRE_NAME
  fire_year <- mtbs_fire$YEAR_MTBS
  fire_size <- mtbs_fire$ACRES
  
  mtbs_fire_ig <- filter(mtbs_point, FIRE_ID == fire_id)
  
  MTT <- raster(here("data", "intermediate", "FB", "TestMTT", "MTT_Inputs", 
                     "MTT_Output", paste0("fire", n, "_ArrivalTime.asc")))
  
  MTT_INT <- raster(here("data", "intermediate", "FB", "TestMTT", "MTT_Inputs", 
                         "MTT_Output", paste0("fire", n, "_INTENSITY.asc")))
  
  #### Get ignition point from centroid of first day
  
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
  
  mtbsBS_year <- raster(here("data", "raw", "MTBS", "MTBS_BSmosaics", fire_year, paste0("mtbs_CONUS_", fire_year, ".tif")))
  
  # intersect MTBS w/FACTS
  mtbs_facts_comp_int <- st_intersection(mtbs_fire, facts_comp)
  mtbs_facts_incomp_int <- st_intersection(mtbs_fire, facts_incomp)
  
  mtbs_facts_comp_int_filt <- as.data.frame(subset(mtbs_facts_comp_int, YEAR < YEAR_MTBS | YEAR == YEAR_MTBS & MONTH < MONTH_MTBS)) %>%
    filter(YEAR > 2004) %>%
    filter(YEAR_MTBS - YEAR <= 10)
  
  mtbs_facts_incomp_int_filt <- as.data.frame(subset(mtbs_facts_incomp_int, YEAR_PLANNED < YEAR_MTBS | YEAR_PLANNED == YEAR_MTBS & MONTH_PLANNED < MONTH_MTBS)) %>%
    filter(YEAR_PLANNED > 2004)
  
  mtbs_facts_int_filt <- rbind(mtbs_facts_comp_int_filt, mtbs_facts_incomp_int_filt)
  
  treatment_IDs <- unique(mtbs_facts_int_filt$ROW_ID)
  
  fire_extent_poly <- st_as_sfc(st_bbox(fire))
  
  # Find intersecting fuel treatments
  intersecting_facts <- filter(facts, ROW_ID %in% treatment_IDs)
  intersecting_facts <- filter(intersecting_facts, COMPLETE == 0) # Only look at incomplete projects for now
  
  if (nrow(intersecting_facts) == 0){
    next
  }
  
  # Get the boundary of the fire perimeter
  fire_boundary <- st_boundary(mtbs_fire)
  
  # Convert fire boundary to multipoint for easier distance calculations
  boundary_points <- st_cast(fire_boundary, "POINT")
  
  # Calculate the distances from the ignition point to each boundary point
  distances <- st_distance(mtbs_fire_ig, boundary_points)
  
  # Find the maximum distance
  max_distance <- max(as.numeric(distances))  # Extract numeric value from the distance object
  
  # add 3 kilometers to max distance
  
  max_distance <- max_distance + 3000
  
  # Define the number of rays
  num_rays <- 24
  
  # Create a sequence of angles (in radians)
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
  distance_bins <- seq(0, max_distance, by = 500)  # 500 meters = 0.5 km
  
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
  
  
  #### Get the dominant angle of spread for each direction
  
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
  
  
  # Now you have rays divided by direction and distance bins
  # You can save this as a shapefile or use it in further analysis
  
  Grids <- direction_distance_sf
  Grids$Grid_ID <- 1:nrow(Grids)
  Grids_Cent <- st_centroid(Grids)
  
  #### Intersect with National Forests
  
  Grids_NF_int <- st_intersection(Grids_Cent, USFS_NF)
  
  grids_usfs_nf_int <- as.data.frame(Grids_NF_int) %>%
    group_by(Grid_ID) %>%
    summarise(USFS_NF = 1) %>%
    dplyr::select(Grid_ID, USFS_NF)
  
  Grids <- merge(Grids, grids_usfs_nf_int, by = "Grid_ID", all.x = T)
  Grids[is.na(Grids$USFS_NF), "USFS_NF"] <- 0
  
  #### Intersect with Wilderness Areas
  
  Grids_NF_int <- st_intersection(Grids_Cent, WAs)
  
  grids_usfs_wa_int <- as.data.frame(Grids_NF_int) %>%
    group_by(Grid_ID) %>%
    summarise(Wilderness = 1) %>%
    dplyr::select(Grid_ID, Wilderness)
  
  Grids <- merge(Grids, grids_usfs_wa_int, by = "Grid_ID", all.x = T)
  Grids[is.na(Grids$Wilderness), "Wilderness"] <- 0
  
  # Calculate Acres in Grid
  Grids$Grid_Acres <- as.numeric(st_area(Grids)/4046.86)
  
  # Group all FACTS treatments by Activity Unit
  
  intersecting_facts_grouped <- intersecting_facts %>%
    group_by(ACTIVITY_UNIT_CN) %>%
    summarise(geometry = st_union(geometry), 
              TREAT_CAT = case_when(length(unique(TREATMENT_CAT)) == 1 ~ first(TREATMENT_CAT),
                                    length(unique(TREATMENT_CAT)) > 1 ~ "Rx & Mech"), 
              YEAR_COMPLETE = max(YEAR, na.rm = T), 
              TREAT_TYPES = paste(unique(TREATMENT_TYPE), collapse = ", "), 
              PLANNED_ACRES = sum(NBR_UNITS_PLANNED, na.rm = TRUE),
              ACCOMPLISHED_ACRES = sum(NBR_UNITS_ACCOMPLISHED, na.rm = TRUE),
              COMPLETE = ifelse(ACCOMPLISHED_ACRES / PLANNED_ACRES < 0.10, 0, 1))
  
  act_treat_complete <- dplyr::select(as.data.frame(intersecting_facts_grouped), ACTIVITY_UNIT_CN, COMPLETE) %>% dplyr::rename(TREAT_ID = ACTIVITY_UNIT_CN)
  
  # Calculate acres treated for each treatment
  intersecting_facts_grouped$Acres_Treated_Treatment <- as.numeric(st_area(intersecting_facts_grouped)/4046.86)
  
  Grids_Treated_Int <- st_intersection(Grids, intersecting_facts_grouped)
  Grids_Treated_Int$Acres_Treated_Int <- as.numeric(st_area(Grids_Treated_Int)/4046.86)
  
  Grids_Treated <- Grids_Treated_Int %>%
    as.data.frame() %>%
    group_by(Grid_ID) %>%
    summarise(Acres_Treated = sum(Acres_Treated_Int), 
              Grid_Acres = dplyr::first(Grid_Acres), 
              TREAT_CAT = case_when(length(unique(TREAT_CAT)) == 1 ~ first(TREAT_CAT),
                                    length(unique(TREAT_CAT)) > 1 ~ "Rx & Mech"), 
              TREAT_SIZE = sum(Acres_Treated_Treatment), 
              YEAR_COMPLETE = max(YEAR_COMPLETE, na.rm = T), 
              TREAT_ID  = dplyr::first(ACTIVITY_UNIT_CN), 
              TREAT_TYPES = paste(unique(TREAT_TYPES), collapse = ", "))  %>% # 
    mutate(Pct_Treated = Acres_Treated/Grid_Acres,
           Treated = ifelse(Pct_Treated > 0.50, 1, 0)) %>%
    dplyr::select(Grid_ID, Pct_Treated, Treated, TREAT_CAT, TREAT_SIZE, YEAR_COMPLETE, TREAT_ID, TREAT_TYPES)
  
  Grids_Treated <- merge(Grids_Treated, act_treat_complete, by = "TREAT_ID", all.x = T)
  
  Grids <- merge(Grids, Grids_Treated, by = "Grid_ID", all.x = T)
  
  
  Grids[is.na(Grids$Treated), "Treated"] <- 0
  Grids[is.na(Grids$Pct_Treated), "Pct_Treated"] <- 0
  
  Grids_lag_treated <- as.data.frame(Grids) %>%
    group_by(direction) %>%   # Group by direction (fire identifier)
    arrange(distance_bin_end) %>%         # Ensure the data is sorted by distance bin
    mutate(Treated_LAG = lag(Treated, n = 1, default = 0)) %>%  # Create the lagged burned variable
    ungroup() %>%
    dplyr::select(direction, distance_bin_end, Treated_LAG) %>%
    arrange(direction)
  
  Grids <- merge(Grids, Grids_lag_treated, 
                 by = c("direction", "distance_bin_end"), all.x = T)
  
  
  
  
  Grids <- mutate(Grids, distance_bin_end = distance_bin_end/500)
  
  ## Create indicators for whether a direction is a treated direction or not. Get distance bin of closest treatment in a direction - L_TAU
  
  Direction_L_Taus <- as.data.frame(Grids) %>%
    group_by(direction) %>%
    filter(Treated == 1) %>%
    summarise(L_TAU = min(as.numeric(distance_bin_end)),
              Treated_Dir = 1)
  
  Grids <- merge(Grids, Direction_L_Taus, by = "direction", all.x = T)
  Grids <- mutate(Grids, AFT_TREAT = ifelse(distance_bin_end > L_TAU & Treated == 0, 1, 0))
  
  Grids[is.na(Grids$Treated_Dir), "Treated_Dir"] <- 0
  
  Grids <- mutate(Grids, PRE = ifelse(distance_bin_end < L_TAU, 1, 0),
                  POST = ifelse(distance_bin_end >= L_TAU, 1, 0))
  Grids[is.na(Grids$PRE), "PRE"] <- 0
  Grids[is.na(Grids$POST), "POST"] <- 0
  
  Grids <- mutate(Grids, time_to_treat = ifelse(Treated_Dir == 1, as.numeric(distance_bin_end) - L_TAU, 0))
  
  ## Define Treatment category, size, and time from treatment based on the value for the first treatment to intersect a direction
  
  Treat_Directions <- as.data.frame(Grids) %>%
    filter(Treated_Dir == 1 & time_to_treat == 0) %>%
    group_by(direction) %>%
    summarise(TREAT_CAT = dplyr::first(TREAT_CAT), 
              TREAT_SIZE = dplyr::first(TREAT_SIZE), 
              YEAR_COMPLETE = dplyr::first(YEAR_COMPLETE), 
              TREAT_TYPES = dplyr::first(TREAT_TYPES))
  
  Grids <- subset(Grids, select = -c(TREAT_CAT, TREAT_SIZE, YEAR_COMPLETE, TREAT_TYPES))
  
  Grids <- merge(Grids, Treat_Directions, by = "direction", all.x = T)
  
  #### Get Day of Burn & Burn Severity
  
  Grids_Cent <- st_centroid(Grids)
  
  Grids_Cent_reproj <- st_transform(Grids_Cent, crs = st_crs(fire))
  
  Grids$DAY_BURN <- raster::extract(fire, Grids_Cent_reproj) # Extract the day a grid cell burns
  
  Grids_Fire_Int <- st_intersection(Grids, mtbs_fire)
  Grids_Fire_Int$Fire_Grid_Acres <- as.numeric(st_area(Grids_Fire_Int)/4046.86)
  
  Grids_Fire_Ints <- Grids_Fire_Int %>%
    as.data.frame() %>%
    group_by(Grid_ID) %>%
    summarise(Fire_Grid_Acres = sum(Fire_Grid_Acres), 
              Grid_Acres = dplyr::first(Grid_Acres)) %>%
    mutate(PCT_BURN = Fire_Grid_Acres/Grid_Acres) %>%
    dplyr::select(Grid_ID, PCT_BURN)
  
  Grids <- merge(Grids, Grids_Fire_Ints, by = "Grid_ID", all.x = T)
  Grids[is.na(Grids$PCT_BURN), "PCT_BURN"] <- 0
  
  # Code a grid as burned if the centroid intersects with the fire
  Grids_Cent_Fire_Int <- st_intersection(Grids_Cent, mtbs_fire)
  Grids_Cent_Fire_Int <- as.data.frame(Grids_Cent_Fire_Int) %>%
    group_by(Grid_ID) %>%
    summarise(BURN = 1)
  
  Grids <- merge(Grids, Grids_Cent_Fire_Int, by = "Grid_ID", all.x = T)  
  Grids[is.na(Grids$BURN), "BURN"] <- 0
  
  
  Fire_Direction_Distance_Burned <- as.data.frame(Grids) %>%
    group_by(direction) %>%   # Group by direction (fire identifier)
    arrange(distance_bin_end) %>%         # Ensure the data is sorted by distance bin
    mutate(BURN_LAG = lag(BURN, n = 1, default = 1)) %>%  # Create the lagged burned variable
    ungroup() %>%
    dplyr::select(direction, distance_bin_end, BURN_LAG) %>%
    arrange(direction)
  
  Grids <- merge(Grids, Fire_Direction_Distance_Burned, 
                 by = c("direction", "distance_bin_end"), all.x = T)
  
  Grids$YEAR <- fire_year
  Grids$FIRE_SIZE <- fire_size
  
  # Calcualte Year Since Treatment
  
  Grids <- mutate(Grids, YEARS_SINCE_TREAT = YEAR - YEAR_COMPLETE)
  
  # For all the grids find the raster that it is closest to
  
  Grids_No_Burn <- st_transform(Grids_Cent, crs(fire))
  
  grid_points <- as(Grids_No_Burn, "Spatial")
  
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
  nearest_indices <- st_nn(st_as_sf(Grids_No_Burn), raster_points, returnDist = TRUE)
  
  # Extract the nearest raster values
  nearest_values <- sapply(nearest_indices$nn, function(idx) raster_points$raster_value[idx])
  
  # Add the raster values to the original Grids_Cent data frame
  Grids_No_Burn <- Grids_No_Burn %>%
    mutate(DAY_BURN_CLOSEST = as.numeric(nearest_values))
  
  Grids <- merge(Grids, dplyr::select(as.data.frame(Grids_No_Burn), Grid_ID, DAY_BURN_CLOSEST), by = "Grid_ID", all.x = T)
  
  # Define date burned as the date where the closest grid burns
  Grids <- mutate(Grids, DATE_BURNED = as.Date(paste0(YEAR, "-", DAY_BURN_CLOSEST), format = "%Y-%j")) 
  
  Grids_reproj <- st_transform(Grids, crs = st_crs(mtbsBS_year))
  
  # Extract MTBS Burn Severity for each grid
  Grids$BURN_SEV <- as.numeric(exact_extract(mtbsBS_year, Grids_reproj, "mean", progress = F))
  
  Grids <- mutate(Grids, BURN_SEV0 = ifelse(is.na(BURN_SEV) == T, 0, BURN_SEV)) # BURN_SEV0 codes non-burned grids as zeros
  
  #### Calculate Medium-High Burn Severity
  
  burn_sev_34_pct <- function(df) {
    valid <- !is.na(df$value)
    num <- sum(df$coverage_fraction[valid & df$value %in% c(3, 4)])
    denom <- sum(df$coverage_fraction[valid])
    if (denom == 0) return(NA) else return(num / denom)
  }
  
  
  Grids$BURN_SEV_MED_HIGH_PCT <- exact_extract(
    mtbsBS_year,
    Grids_reproj,
    fun = burn_sev_34_pct,
    summarize_df = TRUE,
    progress = FALSE
  )
  
  Grids[is.na(Grids$BURN_SEV_MED_HIGH_PCT), "BURN_SEV_MED_HIGH_PCT"] <- 0
  
  Grids$FIRE_ID <- fire_id
  Grids$FIRE_NAME <- fire_name
  
  Grids <- Grids[, c(setdiff(names(Grids), "geometry"), "geometry")] # make sure geometry is last column
  Grids <- dplyr::rename(Grids, distance_bin = distance_bin_end)
  
  
  ## Do treated directions ever hit treatment?
  
  Treat_Dir_Hit <- as.data.frame(Grids) %>%
    filter(Treated_Dir == 1) %>%
    group_by(direction) %>%
    summarise(L_TAU = dplyr::first(L_TAU), 
              Max_Dist_Burn = max(distance_bin[BURN == 1])) %>%
    mutate(HIT_TREAT = ifelse(L_TAU <= Max_Dist_Burn + 1, 1, 0)) %>%
    dplyr::select(direction, HIT_TREAT)
  
  Grids <- merge(Grids, Treat_Dir_Hit, by = "direction", all.x = T)
  
  Grids[is.na(Grids$HIT_TREAT), "HIT_TREAT"] <- 0
  
  #### For each treatment-direction get the distance bin it first hits the treatment - "DIST_TREAT"
  
  Treat_Dist_Hit <- as.data.frame(Grids) %>%
    filter(Treated_Dir == 1) %>%
    group_by(direction, TREAT_ID) %>%
    summarise(DIST_TREAT = min(distance_bin[Treated == 1])) %>%
    dplyr::select(direction,TREAT_ID, DIST_TREAT)
  
  Treat_Dist_Hit <- Treat_Dist_Hit[!is.na(Treat_Dist_Hit$TREAT_ID), ]
  
  Grids <- merge(Grids, Treat_Dist_Hit, by = c("direction", "TREAT_ID"), all.x = T)
  
  Grids[is.na(Grids$DIST_TREAT), "DIST_TREAT"] <- 0
  
  
  ## Get the distance bins for each treatment that intersects with a fire (conditional not yet being extinguished - BURN_LAG == 1) in a given direction
  
  Direction_Distance_Treatments <- as.data.frame(Grids) %>%
    filter(!is.na(TREAT_ID)) %>%
    group_by(direction, TREAT_ID) %>%
    filter(Treated == 1) %>%
    summarise(Min_Dist = min(as.numeric(distance_bin)),
              BURN_LAG = max(BURN_LAG), 
              HIT_TREAT = max(HIT_TREAT)) %>%
    filter(BURN_LAG != 0) %>%
    group_by(direction) %>%
    summarise(DIST_HIT_TREATS = paste(unique(Min_Dist), collapse = ", "))
  
  Grids <- merge(Grids, Direction_Distance_Treatments, by = "direction", all.x = T)
  
  #### For each Grid get an indicator if a grid cell has previously burned in the last 10 years
  
  mtbs_previous_10yr <- filter(mtbs, fire_year > YEAR_MTBS) %>% filter(fire_year - YEAR_MTBS <= 10)
  
  
  Grids_Burned_Int <- st_intersection(Grids, mtbs_previous_10yr)
  Grids_Burned_Int$Acres_Burned_Int <- as.numeric(st_area(Grids_Burned_Int)/4046.86)
  
  Grids_Burned <- Grids_Burned_Int %>%
    as.data.frame() %>%
    group_by(Grid_ID) %>%
    summarise(Previous_Acres_Burned = sum(Acres_Burned_Int), 
              Grid_Acres = dplyr::first(Grid_Acres))  %>%
    mutate(Pct_Burned_10Y = Previous_Acres_Burned/Grid_Acres,
           PREV_BURN_10Y = ifelse(Pct_Burned_10Y > 0, 1, 0)) %>%
    dplyr::select(Grid_ID, Pct_Burned_10Y, PREV_BURN_10Y)
  
  Grids <- merge(Grids, Grids_Burned, by = "Grid_ID", all.x = T)
  
  Grids[is.na(Grids$Pct_Burned_10Y), "Pct_Burned_10Y"] <- 0
  Grids[is.na(Grids$PREV_BURN_10Y), "PREV_BURN_10Y"] <- 0
  
  
  #### Create fire-direction, fire-distance, fire-direction-distance FEs, & Burn Right Away Indicator (i.e. if fire intersects with treatment right away)
  
  Grids <- mutate(Grids, direction_distance_fire_FE = as.factor(paste(direction, distance_bin, FIRE_ID, sep = "-")),
                  direction_fire_FE = paste(direction, FIRE_ID, sep = "-"),
                  distance_fire_FE = paste(distance_bin, FIRE_ID, sep = "-"), 
                  Burn_Right_Away = ifelse(is.na(L_TAU) == F & L_TAU == 1, 1, 0))
  
  #### Extract MTT
  
  Grids_reproj <- st_transform(Grids, crs = st_crs(MTT))
  
  Grids$MTT <- as.numeric(exact_extract(MTT, Grids_reproj, "mean", progress = F))
  
  Fire_Direction_Distance_MTT <- as.data.frame(Grids) %>%
    group_by(direction) %>%   # Group by direction (fire identifier)
    arrange(distance_bin) %>%         # Ensure the data is sorted by distance bin
    mutate(MTT_LAG = lag(MTT, n = 1, default = NA)) %>%  # Create the lagged burned variable
    ungroup() %>%
    dplyr::select(direction, distance_bin, MTT_LAG) %>%
    arrange(direction, distance_bin)
  
  Grids <- merge(Grids, Fire_Direction_Distance_MTT, 
                 by = c("direction", "distance_bin"), all.x = T)  
  
  Grids <- mutate(Grids, DELTA_MTT = ifelse(distance_bin == 1, MTT, MTT - MTT_LAG))
  
  ## Delta MTT Missing Indicator - missing if there is no MTT in previous cell do to missing fuel or if arrival time is lower than in focal cell
  
  Grids <- mutate(Grids, NA_DELTA_MTT = ifelse(is.na(DELTA_MTT) == T | DELTA_MTT < 0, 1, 0))
  
  #### Extract Fireline Intensity
  
  Grids_reproj <- st_transform(Grids, crs = st_crs(MTT_INT))
  
  Grids$FIRE_INTENSITY <- as.numeric(exact_extract(MTT_INT, Grids_reproj, "mean", progress = F))
  
  Grids <- mutate(Grids, LOG_FIRE_INTENSITY = log(FIRE_INTENSITY))
  
  #### Create ID
  
  Grids <- mutate(Grids, ID = paste(Grid_ID, FIRE_ID, sep = "-"))
  Grids <- subset(Grids, select = -c(distance_bin_start, Grid_ID, YEAR_COMPLETE))
  
  Grids <- Grids %>%
    dplyr::select(ID, everything())
  
  SpatialOLS_Grids <- rbind(SpatialOLS_Grids, Grids)
  
}

#### Merge smoke data to Grids

SpatialOLS_Grids <- merge(SpatialOLS_Grids, smoke_data, by = "FIRE_ID", all.x = T)

SpatialOLS_Grids <- mutate(SpatialOLS_Grids, FIRE_PM25_PER_ACRE = FIRE_PM25/FIRE_SIZE, 
                           FIRE_CO2_PER_ACRE = FIRE_CO2/FIRE_SIZE)

## Save Spatial OLS grids

# dir.create(here("data", "intermediate", "SpatialOLS_Grids_Shapefile_V2"))
# 
# st_write(SpatialOLS_Grids, here("data", "intermediate", "SpatialOLS_Grids_Shapefile_V2", "SpatialOLS_Grids.shp"), append=FALSE)




# Grids <- mutate(Grids, Treatment_Cat = ifelse(Treated_Dir == 1 & POST == 1, "Treated", "Control"))
# 
# ei_tiles = get_tiles(st_buffer(fire_extent_poly, 4000), provider = "Esri.WorldTopoMap", zoom = 11, crop = TRUE)
# 
# min_val <- minValue(fire)
# max_val <- maxValue(fire)
# 
# breaks <- seq(min_val - .5, max_val + .5, 1)
# 
# labels <- paste0("Day ", 1:length(breaks))
# 
# tmap_mode("plot")
# 
# pdf(here("output", "maps", "Burro_Fire_2017_Baseline.pdf"))
# 
# plot1 <- tm_shape(ei_tiles, is.master = T) +
#   tm_rgb() +
#   tm_shape(fire) +
#   tm_raster(breaks = seq(180.5, 190.5, 1), title = "Day of Burn", labels = labels) +
#   tm_shape(mtbs_fire_ig) +
#   tm_symbols(col = "black", size = 0.01) +
#   tm_shape(intersecting_facts) +
#   tm_polygons(col = "TREATMENT_TYPE", palette = c("steelblue", "black"), lwd = 2, alpha = 0.35) +
#   tm_layout(legend.outside = T)

# print(plot1)
# 
# dev.off()
# 
# pdf(here("output", "maps", "Burro_Fire_2017_SpatialDiD.pdf"))
# 
# 
# Grids <- mutate(Grids, Type = case_when(BURN_LAG == 0 | USFS_NF == 0 ~ "Missing",
#                                         is.na(L_TAU) == T & BURN_LAG == 1 & USFS_NF == 1 ~ "Control",
#                                         POST == 1 & BURN_LAG == 1 & USFS_NF == 1 ~ "Treated",
#                                         POST == 0 & is.na(L_TAU) == F & BURN_LAG == 1 & USFS_NF == 1 ~ "Yet-to-be treated"))
# 
# Grids <- mutate(Grids, Treatment_Direction = case_when(Treated_Dir == 1 & USFS_NF == 1 ~ "Treated Direction",
#                                                        Treated_Dir == 0 & USFS_NF == 1 ~ "Control Direction",
#                                                        USFS_NF == 0 ~ "Missing"))
# 
# # # Plot using tmap
# plot2 <- tm_shape(ei_tiles, is.master = T) +
#   tm_rgb() +
#   tm_shape(fire) +
#   tm_raster(breaks = breaks, labels = labels, legend.show = F) +
#   tm_shape(mtbs_fire_ig) +
#   tm_symbols(col = "black", size = 0.01) +
#   tm_shape(intersecting_facts) +
#   tm_polygons(col = "TREATMENT_TYPE", palette = c("steelblue", "black"), lwd = 2, alpha = 0.35, title = "Treatment Type") +
#   tm_shape(Grids) +
#   tm_polygons(col = "Type",
#               palette = c("Control" = "#e9c46a",           # Soft mustard yellow
#                           "Yet-to-be treated" = "#a3d9a5",  # Light green
#                           "Treated" = "#1b9e77",            # Dark green
#                           "Missing" = "#cccccc"),           # Light gray
#               title = "Treatment Status", alpha = 0.45) +
#   tm_layout(legend.outside = T)
# 
# print(plot2)
# 
# dev.off()


# tm_shape(ei_tiles, is.master = T) +
#   tm_rgb() +
#   tm_shape(fire) +
#   tm_raster(breaks = breaks, title = "Day of Burn", labels = labels) +
#   tm_shape(mtbs_fire_ig) +
#   tm_symbols(col = "black", size = 0.01) +
#   tm_shape(intersecting_facts) +
#   tm_polygons(col = "TREATMENT_TYPE", palette = c("steelblue", "black"), lwd = 2, alpha = 0.35, title = "Treatment Type") +
#   tm_shape(Grids) +
#   tm_polygons(col = "Treated", alpha = 0.2) +
#   tm_layout(legend.outside = T)



# Grids_df <- as.data.frame(Grids)
#
# Treated_Directions <- Grids_df %>%
#   group_by(direction) %>%
#   summarise(MAX_BURN_DIST = max(distance_bin),
#             L_TAU = dplyr::first(L_TAU)) %>%
#   mutate(Event_Treated = ifelse(MAX_BURN_DIST >= L_TAU,  1, 0)) %>%
#   dplyr::select(direction, Event_Treated)
#
# Grids_df <- merge(Grids_df, Treated_Directions, by = "direction", all.x = T)
# Grids_df[is.na(Grids_df$Event_Treated), "Event_Treated"] <- 0
# Grids_df <- mutate(Grids_df, time_to_treat = ifelse(Event_Treated == 1, as.numeric(distance_bin) - L_TAU, 0))
#
#
#
# Grids_df_filt <- filter(Grids_df, FIRE_NAME == "BURRO" & BURN_LAG == 1)
# Grids_df_filt <- mutate(Grids_df_filt, ANY_TREATED = ifelse(time_to_treat >= 0, 1, 0))
# Grids_df_filt <- mutate(Grids_df_filt, Treat_Post = ifelse(Treated_Dir == 1 & POST == 1, 1, 0))
#
# mod_twfe_1 = feols(BURN ~ i(time_to_treat, Event_Treated, ref = -1) |                    ## Other controls
#                      direction + distance_bin,                             ## FEs                         ## Clustered SEs
#                    data = Grids_df_filt)
#
# OLS_Burned2 <- feols(Burned ~ ANY_TREATED, Grids_df_filt, weights = ~Grid_Acres)
#
#
# iplot(mod_twfe_1,
#       xlab = 'Distance to first treatment (0.5 kilometers)',
#       ylab = "Burn Probability",
#       main = 'Event study - Fire-Direction & Fire-Distance FEs')



################ Calculate Grid Characteristics


pacman::p_load(dplyr,sf, tmap, magrittr, rnaturalearth, rnaturalearthdata, ggplot2, maps, lwgeom, rgeos, raster, stars, haven, stargazer, quantmod, lubridate, tidyr, ggpubr,
               rgdal, exactextractr, tictoc, terra, gtools, here, fixest, modelsummary, readr, rdrobust, prism, parallel,tmaptools,
               OpenStreetMap, maptiles, gifski, stringr, tidycensus)


OLS_Grids <- SpatialOLS_Grids

#### Load in USFS Roads & US Highways

USFS_Roads <- st_read(here("data", "raw", "USFS", "Roads", "S_USA.RoadCore_FS.shp")) %>% st_transform(crs = 5070)

US_Highways <- st_read(here("data", "raw", "Highways", "tl_2016_us_primaryroads.shp")) %>% st_transform(crs = 5070)

Western_WUI <- st_read(here("data", "intermediate", "WUI", "western_wui_2010.shp")) %>% st_transform(crs = 5070)

###### Fuel Type & Fire Regime Characteristics (MFRI)

MFRI <- raster(here("data", "raw", "LANDFIRE", "US_140_MFRI", "Tif", "us_140mfri.tif"))

FBFM13 <- raster(here("data", "raw", "LANDFIRE", "US_105_FBFM13", "Tif", "us_105fbfm13.tif")) # Fuel Model Type in LANDFIRE 2001

###### Load LANDFIRE Data

# Load Slope percent, aspect, & Elev

slope <- raster(here("data", "raw", "LANDFIRE", "LF2020_SlpD_220_CONUS", "Tif", "LC20_SlpD_220.tif"))
aspect <- raster(here("data", "raw", "LANDFIRE", "LF2020_Asp_220_CONUS", "Tif", "LC20_Asp_220.tif"))
elev <- raster(here("data", "raw", "LANDFIRE", "LF2020_Elev_220_CONUS", "Tif", "LC20_Elev_220.tif"))
TRI <- terrain(elev, opt = "TRI")

top_landfire <- stack(slope, aspect, elev, TRI) %>% terra::rast()


# Load Exisiting Vegetation Type

EVT <- raster(here("data", "raw", "LANDFIRE", "LF2020_EVT_220_CONUS", "Tif", "LC20_EVT_220.tif"))

EVT_levels <- as.data.frame(levels(EVT))
unique(EVT_levels$EVT_PHYS)

filter(EVT_levels, EVT_PHYS == "Open Water") # Will use "Open Water" to determine if a grid is inside of a lake this corresponds to a value == 7292

###### Houses & Structures

Buildings <- raster(here("data", "raw", "CommunitiesRisk", "CONUS", "BuildingCount_CONUS.tif"))

HUCount <- raster(here("data", "raw", "CommunitiesRisk", "CONUS", "HUCount_CONUS.tif"))

###### Fire Suppression Effort

QAQC <- st_read(here("data", "raw", "QAQC_lines", "qaqc_lines.shp")) %>% st_transform(crs = 5070)

unique(QAQC$FeatureCat)

## LAT Drops

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

## ACS Data

# ACS_2016_2020 <- read_csv(here("data", "raw", "ACS", "ACS_2016_2020_csv", "nhgis0001_ds249_20205_blck_grp.csv"))
# ACS_2016_2020_sf <- st_read(here("data", "raw", "ACS", "ACS_2016_2020_shape", "US_blck_grp_2020.shp"))
# 
# ACS_2016_2020 <- dplyr::select(ACS_2016_2020, GISJOIN, AMWBE001) %>% dplyr::rename(MED_HOUSE_VAL = AMWBE001)
# ACS_2016_2020_sf <- dplyr::select(ACS_2016_2020_sf, GISJOIN)
# 
# ACS_2016_2020_sf <- merge(ACS_2016_2020_sf, ACS_2016_2020, by = "GISJOIN") %>% st_transform(crs = 5070)



######## Calculate Observable Characteristics of Grids


#### Get Weather Characteristics on day of burn

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

#### ii) Topographic Variables

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

#### iii) Fire Regime Characteristics (MFRI) & Fuel Model Type

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

# OLS_Grids$FBFM13 <- as.numeric(raster::extract(FBFM13, OLS_Grids_Centroids))

rm(FBFM13)

# Will categorize each grid into either "Grass", "Shrub", "Timber", "Slash", or "Other" Fuel Types based on - https://www.fs.usda.gov/rm/pubs/rmrs_gtr175/rmrs_gtr175_367_396.pdf

OLS_Grids <- mutate(OLS_Grids,
                    FuelType_2001 = case_when(FBFM13 == 1 | FBFM13 == 2 | FBFM13 == 3 ~ "Grass",
                                              FBFM13 == 4 | FBFM13 == 5 | FBFM13 == 6 | FBFM13 == 7 ~ "Shrub",
                                              FBFM13 == 8 | FBFM13 == 9 | FBFM13 == 10 ~ "Timber",
                                              FBFM13 == 11 | FBFM13 == 12 | FBFM13 == 13 ~  "Slash",
                                              FBFM13 == 91 | FBFM13 == 92 | FBFM13 == 93 | FBFM13 == 98 | FBFM13 == 99 | FBFM13 == -9999 ~ "Other"))


toc()


OLS_Grids_Centroids_reproj <- st_transform(OLS_Grids_Centroids, crs = st_crs(EVT))

tic()

temp <- mclapply(
  1:(num_cores + 1),
  function(i) as.numeric(raster::extract(EVT, filter(OLS_Grids_Centroids_reproj, group_id == i))),
  mc.cores = (num_cores + 1)
)
OLS_Grids$EVT <- unlist(temp)

OLS_Grids <- mutate(OLS_Grids, WATER = ifelse(is.na(EVT) == F & EVT == 7292, 1, 0))

#### iv) Calculate Distance to USFS Road, US Highway, and WUI (or structures later on)

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

OLS_Grids_reproj <- st_transform(OLS_Grids, crs = st_crs(Buildings))

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


# ## For grids with Non-zero buildings or structures calculate the median housing value from ACS
# 
# OLS_Grids_struc <- filter(OLS_Grids, HU_Count > 0 | Struc_Count > 0) %>% st_centroid()
# 
# Grids_ACS_Int <- st_intersection(OLS_Grids_struc, ACS_2016_2020_sf) %>%
#   as.data.frame() %>%
#   dplyr::select(ID, MED_HOUSE_VAL)
# 
# OLS_Grids <- merge(OLS_Grids, Grids_ACS_Int, by = "ID", all.x = T)
# 
# OLS_Grids <- mutate(OLS_Grids, TOT_HOUSE_VAL = MED_HOUSE_VAL*HU_Count, TOT_STRUC_VAL = MED_HOUSE_VAL*Struc_Count)


#### vi) Fire Suppression Effort Variables

## Get distance to any suppression effort line

temp <- mclapply(
  1:(num_cores + 1),
  function(i) calculate_distance(filter(OLS_Grids, group_id == i), QAQC),
  mc.cores = (num_cores + 1)
)
OLS_Grids$SUP_LINE_DIST <- unlist(temp)

OLS_Grids <- mutate(OLS_Grids, SUP_LINE = ifelse(SUP_LINE_DIST == 0, 1, 0))


## Calculate the length of suppression lines in each grid

# Find the intersections of lines with each grid cell
lines_in_grids <- st_intersection(QAQC, OLS_Grids)

# Calculate the length of each intersected line segment
lines_in_grids <- lines_in_grids %>%
  mutate(segment_length = st_length(geometry))

# Sum the lengths for each grid cell
grid_line_lengths <- as.data.frame(lines_in_grids) %>%
  group_by(ID) %>%
  summarise(SUP_LINE_LEN = as.numeric(sum(segment_length)))

OLS_Grids <- merge(OLS_Grids, grid_line_lengths, by = "ID", all.x = T)

OLS_Grids[is.na(OLS_Grids$SUP_LINE_LEN), "SUP_LINE_LEN"] <- 0

OLS_Grids <- mutate(OLS_Grids, SUP_LINE_INT = SUP_LINE_LEN/Grid_Acres)

## Get Distance to LAT Drops

OLS_Grids_final <- OLS_Grids[0,]

for (yr in seq(2017, 2023, 1)){
  
  print(yr)
  
  # if (yr >= 2020){
  #   next
  # }
  
  OLS_Grids_yr <- filter(OLS_Grids, YEAR == yr)
  OLS_Grids_yr$DIST_LAT <- calculate_distance(OLS_Grids_yr, filter(LAT, YEAR == yr))
  OLS_Grids_final <- rbind(OLS_Grids_final, OLS_Grids_yr)
  
}

OLS_Grids_final <- mutate(OLS_Grids_final, LAT = ifelse(DIST_LAT == 0 & is.na(DIST_LAT) == F, 1, 0))

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


######## Save grids as .csv

OLS_Grids_df <- as.data.frame(OLS_Grids_final) %>% subset(select = -c(geometry, grid_id, group_id))

write_csv(OLS_Grids_df, here("data", "intermediate", "SpatialDiD_Grids_v3_24L_Incomp.csv"))