##########################################   01_create_maps.R  ##########################################

################# Purpose: Creates maps of fires used to create Figure 2 and Figure S2.
################# Outputs: i) Figure 2 saved as "Figure2.pdf" stored in "output/figures".
#################         ii) Figure S2 saved as "FigureS2.pdf" in "output/figures".

################# Estimated run time: ~ 9 min

rm(list=ls())

# if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr,sf, tmap, magrittr, rnaturalearth, rnaturalearthdata, ggplot2, maps, lwgeom, raster, stars, haven, stargazer, quantmod, lubridate, tidyr, ggpubr, 
               exactextractr, tictoc, terra, gtools, here, fixest, modelsummary, readr, rdrobust, prism, parallel,tmaptools, 
               OpenStreetMap, maptiles, gifski, purrr, nngeo, stringr, magick, geosphere)


# Set Path
here::i_am("code/08_maps/01_create_maps.R")

# Load Functions 

source(here("code", "functions","tidy_facts.R"))
source(here("code", "functions","calculate_distance.R"))
st_erase = function(x, y) st_difference(x, st_make_valid(st_union(st_combine(y)))) # Taken from here https://r-spatial.github.io/sf/reference/geos_binary_ops.html

tic()

########    Load in Datasets    ########

#### FACTS - Forest Service Fuel Treatments

facts <- st_read(here("data", "raw", "FACTS", "S_USA.Activity_HazFuelTrt_PL.shp"))
facts <- tidy_facts(facts, inf_yr = 2020, crs = 5070)
facts$ROW_ID <- 1:nrow(facts)

# distinguish completed & incomplete projects
facts_comp <- filter(facts, COMPLETE == 1)
facts_incomp <- filter(facts, COMPLETE == 0)

#### Load in Parks Daily Fire Perimeters

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

mtbs_2006_2023 <- filter(mtbs, YEAR_MTBS >= 2006 & YEAR_MTBS <= 2023)

mtbs_2006_2023_df <- as.data.frame(mtbs_2006_2023) %>% dplyr::select(-geometry)

write_csv(tibble(FIRE_ID = mtbs_2006_2023_df$FIRE_ID), here("data", "intermediate", "mtbs_fire_ids_2006_2023.csv"))


#### Firelines

QAQC <- st_read(here("data", "raw", "QAQC", "FTE_containment.shp")) %>% st_transform(crs = 5070)

QAQC$year <- substr(QAQC$incnum1, 1, 4)

as.data.frame(QAQC) %>%
  group_by(year) %>%
  summarise(N_count = n())

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



fires_int_df <- read_csv(here("data", "intermediate", "FACTS_Parks_Fire_List.csv"))


##################     Build Figure 2       ##################   

#### This uses code from "05_panel/01_build_plot_panel.R" to create the plots used in the figure

# i = 1 Burro Fire

###### Code from "01_build_plot_panel.R"

i = 1

n = i

fire_id <- as.character(fires_int_df[i, 1])
ras_id <- as.numeric(fires_int_df[i, 2])

fire <- raster_list[[ras_id]]
mtbs_fire <- filter(mtbs, FIRE_ID == fire_id)
fire_name <- mtbs_fire$FIRE_NAME
fire_year <- mtbs_fire$YEAR_MTBS
fire_size <- mtbs_fire$ACRES

mtbs_fire_ig <- filter(mtbs_point, FIRE_ID == fire_id)

MTT <- raster(here("data", "raw", "FB", "TestMTT", "MTT_Inputs", 
                   "MTT_Output", paste0("fire", n, "_ArrivalTime.asc")))

MTT_INT <- raster(here("data", "raw", "FB", "TestMTT", "MTT_Inputs", 
                       "MTT_Output", paste0("fire", n, "_INTENSITY.asc")))

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
intersecting_facts <- filter(intersecting_facts, COMPLETE == 1) # Only look at complete projects for now

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
            COMPLETE = ifelse(ACCOMPLISHED_ACRES / PLANNED_ACRES < 0.25, 0, 1))

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

###### End of code from "01_build_plot_panel.R"

Grids <- mutate(Grids, Treatment_Cat = ifelse(Treated_Dir == 1 & POST == 1, "Treated", "Control"))

# get the background for the figure

ei_tiles = get_tiles(st_buffer(fire_extent_poly, 4000), provider = "Esri.WorldTopoMap", zoom = 11, crop = TRUE)

# min and max day burned
min_val <- minValue(fire)
max_val <- maxValue(fire)

breaks <- seq(min_val - .5, max_val + .5, 1)

labels <- paste0("Day ", 1:length(breaks))

intersecting_facts_1 <- sf::st_as_sf(intersecting_facts) %>% 
  dplyr::filter(!is.na(TREATMENT_TYPE)) %>%
  dplyr::mutate(
    TREATMENT_TYPE = factor(as.character(TREATMENT_TYPE),
                            levels = c("Broadcast Burn", "Biomass Removal"))
  )

tmap_mode("plot")

pdf(here("output", "maps", "Burro_Fire_2017_Baseline.pdf"))

plot1 <- tm_shape(ei_tiles, is.master = TRUE) +
  tm_rgb() +
  tm_shape(fire) +
  tm_raster(breaks = breaks, title = "Day of burn",
            labels = labels, palette = "YlOrRd") +
  tm_shape(mtbs_fire_ig) +
  tm_symbols(col = "black", size = 0.01) +
  # tm_shape(intersecting_facts_1) +
  # tm_polygons(col = "TREATMENT_TYPE",
  #             palette = c("steelblue", "black"),
  #             lwd = 2, alpha = 0.35,
  #             title = "Treatment type") +
  tm_shape(intersecting_facts_1) +
  tm_polygons(
    fill = "TREATMENT_TYPE",
    fill.scale  = tm_scale(values = c("steelblue", "black")),
    fill.legend = tm_legend(title = "Treatment type"),
    lwd = 2,
    fill_alpha = 0.35
  ) + 
  tm_shape(LAT_2017) + 
  tm_lines(col = "purple", alpha = 0.5,
           legend.title  = "Suppression effort",
           legend.labels = "Suppression effort") + 
  tm_layout(
    legend.outside        = TRUE,
    legend.outside.position = "right",   # or "bottom" if you prefer
    fontfamily            = "Helvetica",
    title.size            = 0.8,
    legend.title.size     = 0.6,
    legend.text.size      = 0.5,
    legend.spacing = 0.2,
    attr.outside          = TRUE
  ) +
  tm_add_legend(type = "line", col = "purple",
                title = "Suppression effort")

print(plot1)

dev.off()


###### Aside - Create .gif for this map for presentations #######

# Define burn days
unique_days <- sort(unique(values(fire)))
unique_days <- unique_days[!is.na(unique_days)]  # remove NA if present

# Create a temporary folder to store images
dir.create(here("output", "frames"))

# Loop through days and save plots
for (i in seq_along(unique_days)) {
  day_cutoff <- unique_days[i]
  
  fire_terra <- rast(fire)  # convert to SpatRaster
  reclass_matrix <- matrix(c(-Inf, day_cutoff, 1,  
                             day_cutoff, Inf, NA), 
                           ncol = 3, byrow = TRUE)
  
  fire_masked <- terra::classify(fire_terra, reclass_matrix)
  fire_progression <- mask(fire_terra, fire_masked)
  
  # Create the map
  map <- tm_shape(ei_tiles, is.master = TRUE) +
    tm_rgb() +
    tm_shape(fire_progression) +
    tm_raster(breaks = breaks, title = "Day of Burn", labels = labels, palette = "YlOrRd") +
    tm_shape(intersecting_facts) +
    tm_polygons(col = "TREATMENT_TYPE", palette = c("steelblue", "black"), lwd = 2, alpha = 0.35, title = "Treatment Type") +
    tm_layout(legend.outside = TRUE) + 
    tm_add_legend(type = "line", col = "purple", title = "Suppression Effort")
  
  # Save PNG for each frame
  tmap_save(map, filename = here("output", sprintf("frames/frame_%02d.png", i)), width = 1200, height = 1000)
}

# Read images and create GIF
frames <- list.files(here("output", "frames"), pattern = "frame_.*\\.png", full.names = TRUE)
img_list <- image_read(frames)
img_gif <- image_animate(image_join(img_list), fps = 1)  # adjust fps as desired
image_write(img_gif, here("output", "Burro_Progression.gif"))

###### End of aside  #######

###### Create Spatial DiD figure 

Grids <- mutate(Grids, Type = case_when(BURN_LAG == 0 ~ "Excluded",
                                        is.na(L_TAU) == T & BURN_LAG == 1 ~ "Control",
                                        POST == 1 & BURN_LAG == 1 ~ "Treated",
                                        POST == 0 & is.na(L_TAU) == F & BURN_LAG == 1 ~ "Yet-to-be treated"))

Grids_1 <- sf::st_as_sf(Grids) %>% 
  dplyr::filter(!is.na(Type)) %>%
  dplyr::mutate(
    Type = factor(as.character(Type),
                  levels = c("Control", "Excluded", "Treated",  "Yet-to-be treated"))
  )


Grids <- mutate(Grids, Treatment_Direction = case_when(Treated_Dir == 1 ~ "Treated Direction",
                                                       Treated_Dir == 0 ~ "Control Direction",
                                                       BURN_LAG == 0 ~ "Excluded"))


pdf(here("output", "maps", "Burro_Fire_2017_SpatialDiD.pdf"))

# # Plot using tmap
plot2 <- tm_shape(ei_tiles, is.master = TRUE) +
  tm_rgb() +
  tm_shape(mtbs_fire) +
  tm_polygons(col = "red", alpha = 0.2) +
  tm_shape(mtbs_fire_ig) +
  tm_symbols(col = "black", size = 0.01) +
  # tm_shape(intersecting_facts) +
  # tm_polygons(col = "TREATMENT_TYPE",
  #             palette = c("black", "steelblue"),
  #             lwd = 2, alpha = 0.35,
  #             title = "Treatment type") +
  tm_shape(intersecting_facts_1) +
  tm_polygons(
    fill = "TREATMENT_TYPE",
    fill.scale  = tm_scale(values = c("steelblue", "black")),
    fill.legend = tm_legend(title = "Treatment type"),
    lwd = 2,
    fill_alpha = 0.35
  ) + 
  tm_shape(Grids_1) +
  tm_polygons(
    fill = "Type",
    fill.scale  = tm_scale(values = c("#e9c46a", "#cccccc", "#1b9e77", "#a3d9a5")),
    fill.legend = tm_legend(title = "Treatment status"),
    fill_alpha = 0.65,
    lwd = 0.55
  ) +
  # tm_shape(Grids) +
  # tm_polygons(
  #   col = "Type",
  #   palette = c(
  #     "Control"          = "#e9c46a",
  #     "Yet-to-be treated"= "#a3d9a5",
  #     "Treated"          = "#1b9e77",
  #     "Excluded"         = "#cccccc"
  #   ),
  #   title = "Treatment status",
  #   alpha = 0.65
  # ) +
  tm_layout(
    legend.outside        = TRUE,
    legend.outside.position = "right",
    fontfamily            = "Helvetica",
    title.size            = 0.8,
    legend.title.size     = 0.6,
    legend.text.size      = 0.5,
    legend.spacing = 0.2,
    attr.outside          = TRUE
  )

print(plot2)

dev.off()


# Save maps to images to later resize
tmap_save(plot1, here("output", "maps", "Burro_Fire_2017_Baseline.png"),
          width = 8, height = 6, units = "in", dpi = 300,
          outer.margins = 0)

tmap_save(plot2, here("output", "maps", "Burro_Fire_2017_SpatialDiD.png"),
          width = 8, height = 6, units = "in", dpi = 300,
          outer.margins = 0)

img_path1 <- here("output", "maps", "Burro_Fire_2017_Baseline.png")
img_path2 <- here("output", "maps", "Burro_Fire_2017_SpatialDiD.png")

# ---- Trim AND match dimensions ----
trim_and_match <- function(path1, path2) {
  
  # 1. Read + trim both images
  img1 <- image_read(path1) %>% image_trim()
  img2 <- image_read(path2) %>% image_trim()
  
  # 2. Get image sizes
  info1 <- image_info(img1)
  info2 <- image_info(img2)
  
  # 3. Choose common dimensions (max width, max height)
  common_width  <- max(info1$width,  info2$width)
  common_height <- max(info1$height, info2$height)
  
  geom <- sprintf("%dx%d", common_width, common_height)
  
  # 4. Extend both images to same geometry, centered
  img1_fixed <- image_extent(img1, geom,  gravity = "center", color = "white")
  img2_fixed <- image_extent(img2, geom,  gravity = "center", color = "white")
  
  # 5. Write them back to disk
  image_write(img1_fixed, path1)
  image_write(img2_fixed, path2)
}

# Run the function on your two images
trim_and_match(img_path1, img_path2)

# Read images 
img1 <- ggdraw() + draw_image(here("output", "maps", "Burro_Fire_2017_Baseline.png")) 
img2 <- ggdraw() + draw_image(here("output", "maps", "Burro_Fire_2017_SpatialDiD.png"))

## Create layout of Figure 2

final_plot <- (img1 | img2) +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag = element_text(
      family = "Helvetica",
      face   = "bold",
      size   = 10
    ),
    panel.spacing = unit(3, "lines")
  )

## Save Figure 2

ggsave(here("output", "figures", "Figure2.pdf"), final_plot,        
       width  = 7.24,          # Science 3-column width
       height = 4.5,    # ~4.5 in
       units  = "in",
       dpi    = 300)


# #### Map of Full Sample
# 
# ## get full set of fires
# 
# MTBS_IDs <- unique(fires_int_df$ID)
# 
# mtbs_sample <- filter(mtbs, FIRE_ID %in% MTBS_IDs) %>% st_transform(crs = 3857)
# 
# usa <- rnaturalearth::ne_states(country = "United States of America", returnclass = "sf") %>%
#   st_transform(crs = 3857)
# 
# Western_States <- c("Washington", "Oregon", "California", "Idaho", "Utah", "Arizona",
#                     "New Mexico", "Montana", "Colorado", "Wyoming", "Nevada")
# 
# usa_western <- filter(usa, woe_name %in% Western_States)
# 
# ei_tiles <- get_tiles(st_bbox(usa_western), provider = "Esri.WorldTopoMap", zoom = 6, crop = TRUE)
# 
# # High-quality map
# tmap_mode("plot")  # Ensures static map for journal-quality output
# 
# # Expand the bounding box by adding a buffer (e.g., 0.05 units)
# burro_buffer <- st_buffer(mtbs_fire_ig, dist = 50000) %>% st_transform(crs = 3857)
# 
# pdf(here("output", "maps", "SampleFires.pdf"))
# 
# plot3 <- tm_shape(usa_western) +
#   tm_polygons(col = "grey90", border.col = "black", lwd = 1) +  # Light grey background, black borders
#   tm_shape(mtbs_sample) +
#   tm_polygons(col = "red3", alpha = 0.6, border.col = "black", lwd = 0.5) +  # Fires in deep red
#   # tm_shape(burro_buffer) +
#   # tm_borders(lwd = 3, col = "steelblue", alpha = 0.8) +
#   tm_layout(title = "",
#             title.position = c("center", "top"),
#             legend.position = c("right", "bottom"),
#             legend.text.size = 1.2,
#             legend.title.size = 1.5,
#             frame = FALSE)  # Removes default frame
# 
# print(plot3)
# 
# dev.off()
# 
# 
# tmap_save(plot3, here("output", "maps", "SampleFires.png"))
# 
# 
# 
# img3 <- ggdraw() + draw_image(here("output", "maps", "SampleFires.png"))
# 
# 
# top_row <- img3 + 
#   plot_annotation(tag_levels = 'a') & 
#   theme(plot.tag.position = c(0.3, 0.98))  # only affects this row
# 
# bottom_row <- img1 | img2 + 
#   plot_annotation(tag_levels = 'a') 
# 
# final_layout <- top_row / bottom_row +
#   plot_layout(heights = c(1.3, 1.1)) +
#   plot_annotation(tag_levels = 'a')
# 
# final_layout 


##################     Build Figure S2       ##################   

#### Examples of different fires & fuel treatments

fire_examples <- c(4, 8, 105, 46)

# i = 4 - Pinal Fire Very Effective, i = 8 Cellar Fire - Effective, i = 105 Cougar Creek - Ineffective, i = 46 - Blow Right Through

plot_list <- list()

for (i in fire_examples){
  
  #### Again using code from "01_build_plot_panel.R"
  
  n = i
  
  if (i == 4){
    title = "Very Effective - 2017 Pinal Fire"
  }
  
  if (i == 8){
    title = "Effective - 2019 Cellar Fire"
  }
  
  if (i == 105){
    title = "Mixed Effectiveness - 2018 Cougar Creek Fire"
  }
  
  if (i == 46){
    title = "Ineffective - 2018 416 Fire"
  }
  
  fire_id <- as.character(fires_int_df[i, 1])
  ras_id <- as.numeric(fires_int_df[i, 2])
  
  fire <- raster_list[[ras_id]]
  mtbs_fire <- filter(mtbs, FIRE_ID == fire_id)
  fire_name <- mtbs_fire$FIRE_NAME
  fire_year <- mtbs_fire$YEAR_MTBS
  fire_size <- mtbs_fire$ACRES
  
  mtbs_fire_ig <- filter(mtbs_point, FIRE_ID == fire_id)
  
  MTT <- raster(here("data", "raw", "FB", "TestMTT", "MTT_Inputs", 
                     "MTT_Output", paste0("fire", n, "_ArrivalTime.asc")))
  
  MTT_INT <- raster(here("data", "raw", "FB", "TestMTT", "MTT_Inputs", 
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
  intersecting_facts <- filter(intersecting_facts, COMPLETE == 1) # Only look at complete projects for now
  
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
              COMPLETE = ifelse(ACCOMPLISHED_ACRES / PLANNED_ACRES < 0.25, 0, 1))
  
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
  
  #### end of code from "01_build_plot_panel.R"
  
  
  ei_tiles = get_tiles(st_buffer(fire_extent_poly, 4000), provider = "Esri.WorldTopoMap", zoom = 11, crop = TRUE)
  
  min_val <- minValue(fire)
  max_val <- maxValue(fire)
  
  breaks <- seq(min_val - .5, max_val + .5, 1)
  
  labels <- paste0("Day ", 1:length(breaks))
  
  LAT_yr <- filter(LAT, YEAR == fire_year)
  
  intersecting_facts <- mutate(intersecting_facts, 
                               TREAT_CAT = case_when(
                                 TREATMENT_CAT == "Fire" ~ "Prescribed Fire", 
                                 TREATMENT_CAT == "Mechanical" ~ "Mechanical")) 
  
  # Ensure "Prescribed Fire" is always first (blue) and "Mechanical" second (green)
  intersecting_facts <- intersecting_facts %>%
    mutate(TREAT_CAT = factor(TREAT_CAT, levels = c("Prescribed Fire", "Mechanical")))
  
  
  # pdf(here("output", "maps", paste0(fire_name, "_Fire_", fire_year, ".pdf")))
  
  dob_breaks <- c(-0.5, 2.5, 5.5, 10.5, 15.5, 25.5, Inf)
  dob_labels <- c("0–1", "2–4", "5–9", "10–14", "15–24", "25+")
  
  min_day <- raster::minValue(fire)
  
  fire_since_ig <- raster::calc(fire, fun = function(x) x - min_day + 1)
  
  intersecting_facts <- sf::st_as_sf(intersecting_facts) %>% 
    dplyr::filter(!is.na(TREAT_CAT)) %>%
    dplyr::mutate(
      TREAT_CAT = factor(as.character(TREAT_CAT),
                         levels = c("Prescribed Fire", "Mechanical"))
    )
  
  fire_buff <- st_buffer(st_transform(mtbs_fire, st_crs(LAT_yr)), dist = 4000) # 10 km
  LAT_clip <- st_intersection(LAT_yr, fire_buff)
  
  plot <- tm_shape(ei_tiles, is.master = T) +
    tm_rgb() +
    tm_shape(fire_since_ig) +
    tm_raster(breaks = dob_breaks, labels = dob_labels,
              palette = "YlOrRd", title = "Days since ignition") +
    tm_shape(mtbs_fire_ig) +
    tm_symbols(col = "black", size = 0.01) +
    tm_shape(intersecting_facts) +
    tm_polygons(
      fill = "TREAT_CAT",
      fill.scale  = tm_scale(values = c("#377EB8", "#41AB5D")),
      fill.legend = tm_legend(title = "Treatment type"),
      lwd = 2,
      fill_alpha = 0.5
    ) +
    tm_shape(LAT_clip) +
    tm_lines(col = "#984EA3", lwd = 1.2) +
    tm_add_legend(type = "line", col = "#984EA3", lwd = 1.2,
                  labels = "Suppression effort", title = "") + 
    tm_layout(
      legend.outside        = TRUE,
      legend.outside.position = "right",
      fontfamily            = "Helvetica",
      # main.title = title,
      # main.title.size = 0.8,
      # main.title.position = "left",
      legend.title.size     = 0.6,
      legend.text.size      = 0.5,
      legend.spacing = 0.2,
      attr.outside          = TRUE
    )
  
  plot_list[[as.character(n)]] <- plot
  
  # print(plot)
  
  # dev.off()
  
  
}

Pinal_Fire_Plot <- plot_list[["4"]]
Cellar_Fire_Plot <- plot_list[["8"]]
Cougar_Creek_Fire_Plot <- plot_list[["105"]]
four16_Fire_Plot <- plot_list[["46"]]



common_inside_bottom <- tm_layout(
  legend.outside = FALSE,
  legend.position = c("left", "bottom"),
  
  # Reserve bottom strip
  inner.margins = c(0.24, 0.02, 0.02, 0.02),
  
  legend.bg.color = "white",
  legend.bg.alpha = 0.95,
  legend.frame = TRUE,
  
  text.fontfamily = "Helvetica",
  legend.title.size = 0.45,
  legend.text.size  = 0.35,
  legend.spacing    = 0.05,
  frame = TRUE
)

lab_style <- function(letter) {
  tm_credits(
    text = letter,
    position = c("left", "top"),
    size = 1.05,
    fontface = "bold",
    col = "black"
  )
}


A <- Pinal_Fire_Plot + lab_style("A") + common_inside_bottom
B <- Cellar_Fire_Plot + lab_style("B") + common_inside_bottom
C <- Cougar_Creek_Fire_Plot + lab_style("C") + common_inside_bottom
D <- four16_Fire_Plot + lab_style("D") + common_inside_bottom

#### Arange four fires using tmap_arrange

fig <- tmap_arrange(A, B, C, D, ncol = 2
                    #,widths = c(0.46, 0.54)  # right column gets more width
)

#### Save Figure S2

tmap_save(
  fig,
  here("output", "figures", "FigureS2.pdf"),
  width = 7.24, height = 7, units = "in",,
  dpi    = 300,
  device = cairo_pdf
)

toc()



