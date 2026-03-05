##########################################    01_build_mtt_landscapes.R  ##########################################

################# Purpose: Creates landscape files, fire ignition shapefiles, input, and command files to run MTT simulations in FlamMap

################# Outputs: i) landscape.tif files for all fires saved in "data/intermediate/MTT_Landscapes" - these will be converted into .lcp files used in FlamMap
#################         ii)  fire ignition shapefiles saved in "data/raw/FB/TestMTT/MTT_Inputs" - these are used in FlamMap simulations
#################         iii) fire .input files saved in "data/raw/FB/TestMTT/MTT_Inputs" - these are used in FlamMap simulations
#################         iv) "simulateCMD.txt" and "simulate.bat" saved in "data/raw/FB/TestMTT/MTT_Inputs" - code to be run in the terminal that runs the MTT simulations

################# Estimated run time: ~45 min

rm(list=ls())

# if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr,sf, tmap, magrittr, rnaturalearth, rnaturalearthdata, ggplot2, maps, lwgeom, raster, stars, haven, stargazer, quantmod, lubridate, tidyr, ggpubr, 
               exactextractr, tictoc, gtools, here, fixest, modelsummary, readr, rdrobust, prism, parallel,tmaptools, 
               OpenStreetMap, maptiles, gifski, purrr, nngeo, stringr, terra)


# Set Path
here::i_am("code/04_mtt/01_build_mtt_landscapes.R")


########    Load in Datasets    ########

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

###### Load LANDFIRE Data

# Load Slope percent, aspect, & Elev

elev <- raster(here("data", "raw", "LANDFIRE", "LF2020_Elev_220_CONUS", "Tif", "LC20_Elev_220.tif"))
slope <- raster(here("data", "raw", "LANDFIRE", "LF2020_SlpD_220_CONUS", "Tif", "LC20_SlpD_220.tif"))
aspect <- raster(here("data", "raw", "LANDFIRE", "LF2020_Asp_220_CONUS", "Tif", "LC20_Asp_220.tif"))

# Fuel Data

FBFM40 <- raster(here("data", "raw", "LANDFIRE", "US_105_FBFM40", "Tif", "us_105fbfm40.tif")) # Fuel Model Type in LANDFIRE 2001
CC <- raster(here("data", "raw", "LANDFIRE", "US_105_CC", "Tif", "us_105cc.tif")) # Fuel Model Type in LANDFIRE 2001
CH <- raster(here("data", "raw", "LANDFIRE", "US_105_CH", "Tif", "us_105ch.tif")) # Fuel Model Type in LANDFIRE 2001
CBH <- raster(here("data", "raw", "LANDFIRE", "US_105_CBH", "Tif", "us_105cbh.tif")) # Fuel Model Type in LANDFIRE 2001
CBD <- raster(here("data", "raw", "LANDFIRE", "US_105_CBD", "Tif", "us_105cbd.tif")) # Fuel Model Type in LANDFIRE 2001

#### CRS used in FlamMap

landscape_example <- raster(here("data", "raw", "landscape_example.tif"))

Flam_CRS <- crs(landscape_example)


###### Create Inputs for MTT Model - loop through all fires

dir.create(here("data", "raw", "FB", "TestMTT", "MTT_Inputs"))
dir.create(here("data", "raw", "FB", "TestMTT", "MTT_Inputs", "MTT_Output"))
dir.create(here("data", "intermediate", "MTT_Landscapes"))


# Initialize an empty vector
skipped_rounds <- c()


######## The goal of this loop is for each fire in our sample to: i) save the first day of ignition burning polygon
########    ii) create the set of plots that define our sample region of interest, iii) within this region
########     crop the relevant landfire variables to create a landscape raster layer that we save.
########    Finally, we iv) calculate the wind speed and direction at the day and location of ignition based on gridMET.


#### This uses code from the function "create_SpatialDiD_Grids" in the script "01_build_plot_panel.R"

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
  
  DOB <- as.Date(mtbs_fire_ig$Ig_Date)
  
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
  
  fire_extent_poly <- st_as_sfc(st_bbox(fire))
  
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
  
  # Now you have rays divided by direction and distance bins
  # You can save this as a shapefile or use it in further analysis
  
  Grids <- direction_distance_sf
  
  # Ensure CRS of the shapefile matches the raster
  Grids <- st_transform(Grids, crs = crs(elev))
  
  # Crop the raster to the extent of the shapefile
  elev_cropped <- crop(elev, st_bbox(Grids)) 
  slope_cropped <- crop(slope, st_bbox(Grids)) 
  aspect_cropped <- crop(aspect, st_bbox(Grids))
  FBFM40_cropped <- crop(FBFM40, st_bbox(Grids))
  CC_cropped <- crop(CC, st_bbox(Grids))
  CH_cropped <- crop(CH, st_bbox(Grids))
  CBH_cropped <- crop(CBH, st_bbox(Grids))
  CBD_cropped <- crop(CBD, st_bbox(Grids))
  
  
  landscape_raster <- stack(elev_cropped, slope_cropped, aspect_cropped, FBFM40_cropped,
                            CC_cropped, CH_cropped, CBH_cropped, CBD_cropped)
  
  #### Reproject tif's and shapefiles to FlamMap crs
  
  landscape_raster_reprojected <- projectRaster(landscape_raster,  crs = Flam_CRS)
  
  # landscape_raster_reprojected_alt <- projectRaster(landscape_raster_alt,  crs = Flam_CRS)
  
  mtbs_fire_ig <- st_transform(mtbs_fire_ig, crs = Flam_CRS)
  mtbs_fire <- st_transform(mtbs_fire, crs = Flam_CRS)
  
  
  
  #### Get Max Wind Speed at Ignition Point on Day of Ignition
  
  vs_yr <- terra::rast(here("data", "raw", "gridMET", paste("vs_", fire_year, ".nc", sep = ""))) # Wind Speed
  th_yr <- terra::rast(here("data", "raw", "gridMET", paste("th_", fire_year, ".nc", sep = ""))) # Wind Direction
  
  Ignition_WindSpeed <- raster::extract(vs_yr, mtbs_fire_ig) %>%
    pivot_longer(cols = -c(ID), names_to = "date_type") %>%
    dplyr::rename(Row_ID = ID) %>%
    mutate(num = as.numeric(str_sub(date_type, -5, -1))) %>%
    mutate(date = num + lubridate::ymd("1900-01-01")) %>%
    filter(date %in% DOB) %>%
    arrange(Row_ID) %>%
    dplyr::rename(WindSpeed = value, DATE_BURNED = date) %>%
    mutate(WindSpeed = 2.23694*WindSpeed) %>% # convert from meters/sec to mph 
    dplyr::select(DATE_BURNED, WindSpeed)
  
  Ignition_WindDirection <- raster::extract(th_yr, mtbs_fire_ig) %>%
    pivot_longer(cols = -c(ID), names_to = "date_type") %>%
    dplyr::rename(Row_ID = ID) %>%
    mutate(num = as.numeric(str_sub(date_type, -5, -1))) %>%
    mutate(date = num + lubridate::ymd("1900-01-01")) %>%
    filter(date %in% DOB) %>%
    arrange(Row_ID) %>%
    dplyr::rename(WindDirection = value, DATE_BURNED = date) %>%
    dplyr::select(WindDirection)
  
  ig_wind_speed_dir <- cbind(Ignition_WindSpeed, Ignition_WindDirection)
  
  ## Save wind speed (WS) and wind direction (WD)
  
  WS <- round(ig_wind_speed_dir$WindSpeed, digits = 0)
  WD <- ig_wind_speed_dir$WindDirection
  
  ## Save Landscape file in "MTT_Landscapes" so they can be converted into a .lcp files
  
  writeRaster(landscape_raster_reprojected, here("data", "intermediate","MTT_Landscapes",  paste0("fire", n, "_landscape.tif")), overwrite = TRUE)
  
  ## Save ignition to MTT_Inputs folder
  
  st_write(ignition_sf, here("data", "raw", "FB", "TestMTT", "MTT_Inputs", paste0("fire",n, "_ig.shp")), append = TRUE)
  
  
  ## Create & write in .input files for each fire
  
  # Define the content of the input file as a character vector
  input_content <- c(
    "ShortTerm-Inputs-File-Version-1",
    "",
    "FUEL_MOISTURES_DATA: 1", 
    "0 6 7 8 60 90", # Standard fuel moisture parameters 
    "",
    paste0("WIND_SPEED: ", WS),
    paste0("WIND_DIRECTION: ", WD),
    "GRIDDED_WINDS_GENERATE: No",
    "GRIDDED_WINDS_RESOLUTION: 200",
    "FOLIAR_MOISTURE_CONTENT: 100",
    "CROWN_FIRE_METHOD: Finney",
    "#CROWN_FIRE_METHOD: ScottReinhardt",
    "NUMBER_PROCESSORS:8",
    "",
    "MTT_RESOLUTION: 150",
    "MTT_SIM_TIME: 0",
    "MTT_TRAVEL_PATH_INTERVAL: 500",
    "MTT_SPOT_PROBABILITY: 0.0",
    "",
    "MTT_FILL_BARRIERS: 0",
    "",
    "NodeSpreadNumLat: 6",
    "NodeSpreadNumVert: 4",
    "#SELECTED MTT OUTPUTS",
    "#MTT_ARRIVAL",
    "#END SELECTED MTT OUTPUTS",
    "#SPOTTING_SEED: -5251"
  )
  
  # Write the content to a file
  writeLines(input_content, here("data", "raw", "FB", "TestMTT", "MTT_Inputs", paste0("fire", n, ".input")))
  
}


#### Now Create "simulateCMD.txt" which is command lines to be run in the TestMTT software.

# Number of fires
N <- nrow(fires_int_df)

# Initialize an empty character vector for the lines
simulate_cmd_lines <- character(N)

# Loop through each fire to create the corresponding line
for (i in 1:N) {
  simulate_cmd_lines[i] <- sprintf(
    "fire%d_landscape.lcp fire%d.input fire%d_ig.shp 0 .\\MTT_Output\\fire%d 1",
    i, i, i, i
  )
}

# Write the lines to "simulateCMD.txt"
writeLines(simulate_cmd_lines, here("data", "raw", "FB", "TestMTT", "MTT_Inputs", "simulateCMD.txt"))


#### Now create "simulate.bat" which is the batch command to run the simulation

# Define the content of the .bat file
bat_content <- "..\\..\\bin\\TestMTT simulateCmd.txt"

# Write the content to "simulate.bat"
writeLines(bat_content, here("data", "raw", "FB", "TestMTT", "MTT_Inputs", "simulate.bat"))