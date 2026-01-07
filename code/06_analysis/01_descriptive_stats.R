##########################################    01_descriptive_stats ##########################################

################# Purpose: Create descriptive statistics Figure 1, impute PM2.5 smoke exposure, and calculate USFS budget imputed footprint cost/acre.

################# Outputs: "SpatialDiD_Grids_L24_K05.csv" - saved in "data/intermediate"
#################           Figure S10 = "Emissions_Exposure_Plot.pdf" saved in "output/figures"
#################           USFS budget imputed footprint cost per acre 
#################           Figure 1 = "Figure1.pdf"

rm(list=ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr,sf, tmap, magrittr, rnaturalearth, rnaturalearthdata, ggplot2, maps, lwgeom, rgeos, raster, stars, haven, stargazer, quantmod, lubridate, tidyr, ggpubr, 
               rgdal, exactextractr, tictoc, terra, gtools, here, fixest, modelsummary, readr, rdrobust, prism, parallel,tmaptools, 
               OpenStreetMap, maptiles, gifski, geosphere, did, cowplot, readxl, scales, gridExtra, patchwork, broom, forcats, magick, stringr)

# Set Path
here::i_am("code/06_analysis/01_descriptive_stats.R")


# Load Functions 

source(here("code", "functions","tidy_facts.R"))
st_erase = function(x, y) st_difference(x, st_make_valid(st_union(st_combine(y)))) # Taken from here https://r-spatial.github.io/sf/reference/geos_binary_ops.html

###### Load in relevant datatset

## Sample of fires

Grids_df <- read_csv(here("data", "temp", "SpatialDiD_Grids_L24_K05.csv"))
MTBS_IDs <- unique(Grids_df$FIRE_ID)


## FRED Elderly Population Data

FRED_65_Over <- read_csv(here("data", "raw", "FRED", "SPPOP65UPTOZSUSA.csv")) %>% 
  dplyr::rename(Pct_Elderly_USA = SPPOP65UPTOZSUSA)

FRED_65_Over <- mutate(FRED_65_Over, YEAR = year(as.Date(observation_date))) %>% dplyr::select(YEAR, Pct_Elderly_USA)

## Load in ICS

ICS_209 <- read.csv(here("data", "raw", "ICS-209", "ics209plus-wildfire", "ics209-plus-wf_incidents_1999to2020.csv"))
ICS_209 <- rename(ICS_209, SUP_COST = PROJECTED_FINAL_IM_COST)

## Adjust Suppression Costs for Inflation
# getSymbols("CPIAUCSL", src='FRED') 

cpi_data <- read_csv(here("data", "raw", "FRED", "CPIAUCSL.csv")) %>%
  mutate(YEAR = year(observation_date))

# Average CPI by year
avg_cpi <- cpi_data %>%
  group_by(YEAR) %>%
  summarise(AVG_CPI = mean(CPIAUCSL, na.rm = TRUE))

# Reference CPI (e.g., 2020)
ref_cpi <- avg_cpi %>%
  filter(YEAR == 2023) %>%
  pull(as.numeric(AVG_CPI))

# Inflation adjustment factor
avg_cpi <- avg_cpi %>%
  mutate(CPI_FACTOR = AVG_CPI / ref_cpi) %>%
  dplyr::rename(START_YEAR = YEAR) %>% 
  dplyr::select(START_YEAR, CPI_FACTOR)

# Join with main data and adjust suppression costs
ICS_209 <- ICS_209 %>%
  left_join(avg_cpi, by = "START_YEAR") %>%
  mutate(SUP_COST_ADJ = SUP_COST / CPI_FACTOR)

# dplyr::select(ICS_209, START_YEAR, SUP_COST, SUP_COST_ADJ)

## Load in USFS NF

USFS_NF <- st_read(here("data", "raw", "USFS", "NationalForests", "S_USA.AdministrativeForest.shp")) %>% st_transform(crs = 5070)

## Group ICS by MTBS Fire ID

ICS_209$FIRE_ID <- gsub("\\(.*\\)", "", ICS_209$LRGST_MTBS_FIRE_INFO)

ICS_209_Grouped <- ICS_209 %>%
  group_by(FIRE_ID) %>%
  summarise(Total_Sup_Cost = sum(SUP_COST_ADJ), 
            Total_Structures_Destroyed = sum(STR_DESTROYED_TOTAL))

ICS_209 %>%
  filter(START_YEAR >= 2017 & START_YEAR <= 2019) %>%  
  summarise(Total_Structures_Destroyed = sum(STR_DESTROYED_TOTAL))

#### FACTS - Forest Service Fuel Treatments

facts <- st_read(here("data", "raw", "FACTS", "S_USA.Activity_HazFuelTrt_PL.shp"))
facts <- tidy_facts(facts, inf_yr = 2023, crs = 5070)
facts$ROW_ID <- 1:nrow(facts)
facts <- mutate(facts, RD_CODE = as.numeric(paste0(ADMIN_REGION_CODE, ADMIN_FOREST_CODE, ADMIN_DISTRICT_CODE)))

states <- c("WA", "OR", "CA", "ID", "NV", "AZ", "MT", "UT", "NM", "WY", "CO")
facts_west <- filter(facts, COMPLETE == 1) %>% filter(STATE_ABBR %in% states)

mean_cost <- mean(facts_west$COST_PER_UOM, na.rm = TRUE)
sd_cost <- sd(facts_west$COST_PER_UOM, na.rm = TRUE)
threshold <- mean_cost + 10*sd_cost

# facts_west_filt <- filter(facts_west, COST_PER_UOM <= threshold)
facts_west_filt <- mutate(facts_west, ADJUSTED_TOTAL_COST = ifelse(COST_PER_UOM <= threshold, 
                                                                   ADJUSTED_TOTAL_COST, 
                                                                   0))
## Load in USA

usa <- rnaturalearth::ne_states(country = "United States of America", returnclass = "sf") %>%
  st_transform(crs = 5070)

# 11 Western States - WA, OR, CA, NV, ID, MT, NM, AZ, UT, CO, WY

usa_west <- filter(usa, name == "Washington" | name == "Oregon" | name == "California" | name == "Nevada" | name == "Idaho" | 
                     name == "Montana" | name == "New Mexico" | name == "Arizona" | name == "Utah" | name == "Colorado" | name == "Wyoming") %>%
  st_transform(crs = 5070)


## Load in MTBS

mtbs <- st_read(here("data", "raw", "MTBS", "MTBS_Burn_Area", "S_USA.MTBS_BURN_AREA_BOUNDARY.shp")) %>%
  st_transform(mtbs, crs = 5070) %>%
  st_make_valid() %>%
  dplyr::rename(YEAR_MTBS = YEAR, MONTH_MTBS = STARTMONTH)

mtbs <- filter(mtbs, FIRE_TYPE == "Unknown" | FIRE_TYPE == "Wildfire") # No Wildland Fire Use Fires


mtbs_point <- st_read(here("data", "raw", "MTBS", "mtbs_fod_pts_data", "mtbs_FODpoints_DD.shp")) %>%
  st_transform(mtbs, crs = 5070) %>%
  st_make_valid()  %>%
  dplyr::rename(FIRE_ID = Event_ID)

#### Load in emissions smoke Data

WFEIS <- read_csv(here("data", "raw", "WFEIS", "WFEIS_data.csv"))


########  Step 1. Impute smoke exposure for fires without PM2.5 exposure estimates

#### What is the correlation between total PM2.5 Emissions from a fire and total population-day PM2.5 exposure?

Fires_PM_Emissions_Exposure <- Grids_df %>%
  group_by(FIRE_ID) %>%
  summarise(FIRE_PM25 = dplyr::first(FIRE_PM25),
            YEAR = dplyr::first(YEAR),
            FIRE_SIZE = dplyr::first(FIRE_SIZE),
            total_pop_smokePM = dplyr::first(total_pop_smokePM))

emission_exposure_plot <- ggplot(
  Fires_PM_Emissions_Exposure,
  aes(x = log(FIRE_PM25), y = log(total_pop_smokePM))
) +
  geom_point(size = 2, alpha = 0.7, color = "gray30") +
  geom_smooth(method = "lm", se = TRUE,
              color = "#0072B2", size = 1.2, fill = "grey80") +
  labs(
    x = expression(bold("Log fire PM"[2.5] * " emissions (" * mu * "g/m"^3 * ")")),
    y = expression(bold("Log total population-day PM"[2.5] * " smoke exposure")),
    title = NULL
  ) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid   = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.title   = element_text(size = 8),   # don’t rely on face="bold" here
    axis.text    = element_text(size = 7),
    axis.ticks   = element_line(linewidth = 0.3),
    plot.title   = element_text(size = 8.5, face = "bold"),
    plot.margin  = margin(3, 3, 3, 3)
  )

final_layout <- emission_exposure_plot +
  plot_layout(guides = "collect") +  # <- collect shared legend
  theme(legend.position = "bottom")  # <- move legend to bottom

ggsave(here("output", "figures", "Emissions_Exposure_Plot.pdf"), final_layout, 
       width  = 7.24,          # Science 3-column width
       height = 4.5,    # ~4.5 in
       units  = "in",
       dpi    = 300)

model <- lm(log(total_pop_smokePM) ~ log(FIRE_PM25), data = Fires_PM_Emissions_Exposure)
summary(model)

# predict for *all* fires (where FIRE_PM25 is observed)
pred_all <- Fires_PM_Emissions_Exposure %>%
  filter(!is.na(FIRE_PM25)) %>%
  mutate(
    predicted_log = predict(model, newdata = .),
    total_pop_smokePM_pred = exp(predicted_log)
  ) %>%
  dplyr::select(FIRE_ID, YEAR, FIRE_SIZE, total_pop_smokePM_pred)

# attach elderly share (by year) and compute predicted deaths/costs/etc
pred_all <- pred_all %>%
  left_join(FRED_65_Over, by = "YEAR") %>%
  mutate(
    Deaths_pred       = (total_pop_smokePM_pred/1000000) * (Pct_Elderly_USA/100) * 0.69,
    EarningsLoss_pred = (total_pop_smokePM_pred/92) * 124.44,
    DeathCost_pred    = Deaths_pred * 11.16 * 1000000,
    
    Deaths_PER_ACRE_pred       = Deaths_pred / FIRE_SIZE,
    DeathCost_PER_ACRE_pred    = DeathCost_pred / FIRE_SIZE,
    EarningsLoss_PER_ACRE_pred = EarningsLoss_pred / FIRE_SIZE
  ) %>%
  dplyr::select(
    FIRE_ID,
    total_pop_smokePM_pred,
    Deaths_pred, EarningsLoss_pred, DeathCost_pred,
    Deaths_PER_ACRE_pred, DeathCost_PER_ACRE_pred, EarningsLoss_PER_ACRE_pred
  )

Grids_df <- Grids_df %>%
  left_join(pred_all, by = "FIRE_ID") %>%
  mutate(
    exposure_predicted = is.na(total_pop_smokePM) & !is.na(total_pop_smokePM_pred),
    
    # fill only if missing
    total_pop_smokePM     = coalesce(total_pop_smokePM,     total_pop_smokePM_pred),
    Deaths                = coalesce(Deaths,                Deaths_pred),
    EarningsLoss          = coalesce(EarningsLoss,          EarningsLoss_pred),
    DeathCost             = coalesce(DeathCost,             DeathCost_pred),
    Deaths_PER_ACRE       = coalesce(Deaths_PER_ACRE,       Deaths_PER_ACRE_pred),
    DeathCost_PER_ACRE    = coalesce(DeathCost_PER_ACRE,    DeathCost_PER_ACRE_pred),
    EarningsLoss_PER_ACRE = coalesce(EarningsLoss_PER_ACRE, EarningsLoss_PER_ACRE_pred)
  )


Grids_df <- arrange(Grids_df, FIRE_ID, direction, distance_bin)

write_csv(Grids_df, here("data", "intermediate", "SpatialDiD_Grids_L24_K05.csv"))




######## Step 2. Create Figure 1 - Map of MTBS fires in sample and time series plots of: 
########                            i) acres burned, ii) suppression costs, iii) acres of fuel treatment, 
########                              and iv) fuel treatment costs in Western U.S. by Western USFS or full Western.


## Merge MTBS with ICS to get costs & structure loss data

mtbs <- merge(mtbs, WFEIS, by = "FIRE_ID", all.x = T)
mtbs <- merge(mtbs, ICS_209_Grouped, by = "FIRE_ID", all.x = T)


# ###### Aside calculate the number of fires intersecting with LAT drops and the number of unique drops ######
# 
# ## LAT Drops
# 
# LAT_2017 <- st_read(here("data", "raw", "LAT", "drops17.shp")) %>% st_transform(crs = 5070)
# LAT_2017$YEAR <- 2017
# LAT_2017$ID <- LAT_2017$Drop_ID1
# LAT_2017 <- dplyr::select(LAT_2017, YEAR, ID)
# 
# LAT_2018 <- st_read(here("data", "raw", "LAT", "drops18.shp")) %>% st_transform(crs = 5070)
# LAT_2018$YEAR <- 2018
# LAT_2018$ID <- LAT_2018$Drop_ID1
# LAT_2018 <- dplyr::select(LAT_2018, YEAR, ID)
# 
# LAT_2019 <- st_read(here("data", "raw", "LAT", "drops19.shp")) %>% st_transform(crs = 5070)
# LAT_2019$YEAR <- 2019
# LAT_2019$ID <- LAT_2019$Drop_ID1
# LAT_2019 <- dplyr::select(LAT_2019, YEAR, ID)
# 
# LAT_2020_2021 <- st_read(here("data", "raw", "LAT", "drops20_21.shp")) %>% st_transform(crs = 5070)
# LAT_2020_2021 <- mutate(LAT_2020_2021, YEAR = localYear1)
# LAT_2020_2021$ID <- LAT_2020_2021$Drop_ID1
# LAT_2020_2021 <- dplyr::select(LAT_2020_2021, YEAR, ID)
# 
# LAT_2022 <- st_read(here("data", "raw", "LAT", "2022_Full_Year_03112024_VLAT.shp")) %>% st_transform(crs = 5070)
# LAT_2022 <- mutate(LAT_2022, YEAR = Year)
# LAT_2022$ID <- LAT_2022$dropID
# LAT_2022 <- dplyr::select(LAT_2022, YEAR, ID)
# 
# LAT_2023 <- st_read(here("data", "raw", "LAT", "2023_Full_Year_03112024_VLAT.shp")) %>% st_transform(crs = 5070)
# LAT_2023 <- mutate(LAT_2023, YEAR = Year)
# LAT_2023$ID <- LAT_2023$dropID
# LAT_2023 <- dplyr::select(LAT_2023, YEAR, ID)
# 
# LAT <- rbind(LAT_2017, LAT_2018, LAT_2019, LAT_2020_2021, LAT_2022, LAT_2023)
# 
# #### Intersect our sample of fires with LAT drops to get the number of WFs and unique drops
# 
# mtbs_sample <- filter(mtbs, FIRE_ID %in% MTBS_IDs)
# 
# lat_summary <- data.frame(
#   Year = integer(),
#   Unique_Fires = integer(),
#   Unique_Drops = integer(),
#   stringsAsFactors = FALSE
# )
# 
# for (yr in seq(2017, 2023, 1)){
#   
#   print(yr)
#   
#   # if (yr >= 2020){
#   #   next
#   # }
#   
#   mtbs_sample_yr <- filter(mtbs_sample, YEAR_MTBS == yr)
#   LAT_yr <- filter(LAT, YEAR == yr)
#   LAT_mtbs_int <- st_intersection(LAT_yr, mtbs_sample_yr)
#   
#   fire_ids <- unique(LAT_mtbs_int$FIRE_ID)
#   drop_ids <- unique(LAT_mtbs_int$ID)
#   
#   lat_summary <- rbind(lat_summary, data.frame(
#     Year = yr,
#     Unique_Fires = length(fire_ids),
#     Unique_Drops = length(drop_ids)
#   ))
#   
# }
# 
# sum(lat_summary$Unique_Fires)
# sum(lat_summary$Unique_Drops)
# 

#### Back to code  ####

## Intersect MTBS with Western U.S. to get list of Western U.S. Fires

ids_usfs <- as.data.frame(st_intersection(mtbs_point, USFS_NF)) %>% group_by(FIRE_ID) %>% summarise(USFS_NF = 1)
ids_western_us <- as.data.frame(st_intersection(mtbs, usa_west)) %>% group_by(FIRE_ID) %>% summarise(WEST = 1)

mtbs <- merge(mtbs, ids_usfs, by = "FIRE_ID", all.x = T)
mtbs <- merge(mtbs, ids_western_us, by = "FIRE_ID", all.x = T)

mtbs[is.na(mtbs$USFS_NF), "USFS_NF"] <- 0
mtbs[is.na(mtbs$WEST), "WEST"] <- 0


mtbs_df <- as.data.frame(mtbs) %>% filter(YEAR_MTBS %in% seq(2006, 2023, 1))

mtbs_df <- mutate(mtbs_df, USFS_WEST = ifelse(USFS_NF == 1 & WEST == 1, 1, 0))

mtbs_region_stats <- mtbs_df %>%
  group_by(YEAR_MTBS) %>%
  summarise(YEARLY_SUP_COST = sum(Total_Sup_Cost, na.rm = T), 
            YEARLY_CO2 = sum(FIRE_CO2, na.rm = T), 
            YEARLY_PM25 = sum(FIRE_PM25, na.rm = T), 
            YEARLY_STRUC = sum(Total_Structures_Destroyed, na.rm = T),
            YEARLY_SUP_COST_USFS = sum(ifelse(USFS_NF == 1, Total_Sup_Cost, NA), na.rm = T),
            YEARLY_CO2_USFS = sum(ifelse(USFS_NF == 1, FIRE_CO2, NA), na.rm = T), 
            YEARLY_PM25_USFS = sum(ifelse(USFS_NF == 1, FIRE_PM25, NA), na.rm = T), 
            YEARLY_STRUC_USFS = sum(ifelse(USFS_NF == 1, Total_Structures_Destroyed, NA), na.rm = T),
            YEARLY_SUP_COST_WEST = sum(ifelse(WEST == 1, Total_Sup_Cost, NA), na.rm = T),
            YEARLY_CO2_WEST = sum(ifelse(WEST == 1, FIRE_CO2, NA), na.rm = T), 
            YEARLY_PM25_WEST = sum(ifelse(WEST == 1, FIRE_PM25, NA), na.rm = T), 
            YEARLY_STRUC_WEST = sum(ifelse(WEST == 1, Total_Structures_Destroyed, NA), na.rm = T),
            YEARLY_SUP_COST_USFS_WEST = sum(ifelse(USFS_WEST == 1, Total_Sup_Cost, NA), na.rm = T),
            YEARLY_CO2_USFS_WEST = sum(ifelse(USFS_WEST == 1, FIRE_CO2, NA), na.rm = T), 
            YEARLY_PM25_USFS_WEST = sum(ifelse(USFS_WEST == 1, FIRE_PM25, NA), na.rm = T), 
            YEARLY_STRUC_USFS_WEST = sum(ifelse(USFS_WEST == 1, Total_Structures_Destroyed, NA), na.rm = T))

mtbs_long2 <- mtbs_df %>%
  filter(WEST == 1) %>%  # Focus only on Western fires
  mutate(
    FIRE_CO2_USFS = ifelse(USFS_NF == 1, FIRE_CO2, 0),
    FIRE_PM25_USFS = ifelse(USFS_NF == 1, FIRE_PM25, 0),
    STRUCT_USFS = ifelse(USFS_NF == 1, Total_Structures_Destroyed, 0),
    SUP_COST_USFS = ifelse(USFS_NF == 1, Total_Sup_Cost, 0),
    ACRES_USFS = ifelse(USFS_NF == 1, ACRES, 0)
  ) %>%
  group_by(YEAR_MTBS) %>%
  summarise(
    CO2_WEST = sum(FIRE_CO2, na.rm = TRUE),
    CO2_USFS = sum(FIRE_CO2_USFS, na.rm = TRUE),
    PM25_WEST = sum(FIRE_PM25, na.rm = TRUE),
    PM25_USFS = sum(FIRE_PM25_USFS, na.rm = TRUE),
    STRUCT_WEST = sum(Total_Structures_Destroyed, na.rm = TRUE),
    STRUCT_USFS = sum(STRUCT_USFS, na.rm = TRUE),
    COST_WEST = sum(Total_Sup_Cost, na.rm = TRUE),
    COST_USFS = sum(SUP_COST_USFS, na.rm = TRUE),
    ACRES_WEST = sum(ACRES, na.rm = TRUE),
    ACRES_USFS = sum(ACRES_USFS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -YEAR_MTBS,
    names_to = c("Metric", "Region"),
    names_sep = "_",
    values_to = "Total"
  ) %>%
  mutate(
    Region = recode(Region, "WEST" = "Western All", "USFS" = "Western USFS"),
    Metric = recode(Metric,
                    "CO2" = "CO₂ Emissions",
                    "PM25" = "PM₂.₅ Emissions",
                    "STRUCT" = "Structures Destroyed",
                    "COST" = "Suppression Costs", 
                    "ACRES" = "Acres Burned")
  )


#### Calculate Acres Burned in USFS

mtbs_USFS_int <- st_intersection(mtbs, USFS_NF)

mtbs_usfs_year_int <- mtbs_USFS_int %>%
  filter(YEAR_MTBS %in% seq(2006,2023, 1)) %>%
  group_by(YEAR_MTBS) %>%
  summarise(geometry = st_union(geometry)) %>%
  arrange(YEAR_MTBS)

mtbs_usfs_year_int$ACRES_BURNED <- as.numeric(st_area(mtbs_usfs_year_int)/4046.86)

mtbs_long2[mtbs_long2$Metric == "Acres Burned" & mtbs_long2$Region == "Western USFS",  "Total"] <- mtbs_usfs_year_int$ACRES_BURNED

## "facts_union_yr" takes the union of all treatments in FACTS to calculate footprint acres treated.
##        to speed up computation it unions treatments at the ranger district level and then aggregates across districts
##        ***Note*** this uses parallel processing 

facts_union_yr <- function(facts_yr){
  
  # Define the function to process each RD_CODE
  process_union <- function(RD) {
    # Filter for the current RD_CODE
    facts_RD <- filter(facts_yr, RD_CODE == RD)
    
    # Check if there are geometries to union
    if (nrow(facts_RD) == 0) {
      return(NULL)  # Skip this RD_CODE if no geometries are found
    }
    
    # Perform the union
    facts_RD_union <- st_union(facts_RD)
    
    # Check if the union resulted in empty geometry
    if (length(facts_RD_union) == 0) {
      return(NULL)  # Skip this RD_CODE if union resulted in empty geometry
    }
    
    # Convert to an sf object and add RD_CODE
    facts_RD_union_sf <- st_sf(geometry = facts_RD_union, RD_CODE = RD)
    
    return(facts_RD_union_sf)
  }
  
  # Use mclapply to process in parallel
  num_cores <- 8  # Use one less core than available
  facts_union_list <- mclapply(RDs, process_union, mc.cores = num_cores)
  
  # Combine all the individual union results into one sf object
  facts_union_yr <- do.call(rbind, facts_union_list)
  
  return(facts_union_yr)
  
}

facts_west <- mutate(facts_west, RD_CODE = as.numeric(paste0(ADMIN_REGION_CODE, ADMIN_FOREST_CODE, ADMIN_DISTRICT_CODE)))
RDs <- unique(facts_west$RD_CODE)

footprint_acres_Rx_Mech <- lapply(seq(2006, 2023), function(i) {
  # Filter facts_filt for the given year i
  
  
  facts_yr <- filter(facts_west, YEAR == i)
  
  # Apply facts_union_yr function to the filtered data
  sum(as.numeric(st_area(facts_union_yr(facts_yr))/4046.86))
})

yearly_acres_treated_USFS <- as.data.frame(facts_west_filt) %>%
  filter(YEAR %in% seq(2006, 2023, 1)) %>%
  group_by(YEAR) %>%
  summarise(YEARLY_COST = sum(ADJUSTED_TOTAL_COST, na.rm = TRUE), 
            ACRES_TREATED = sum(NBR_UNITS_ACCOMPLISHED, na.rm = TRUE))

yearly_acres_treated_USFS$ACRES_TREATED <- as.numeric(unlist(footprint_acres_Rx_Mech)) 


#### Aside Calculate Average Cost/Footprint Acres Treated in whole U.S. from 2011-2020  ####

RDs <- unique(facts$RD_CODE)

facts_2011_2020 <- filter(facts, YEAR %in% seq(2011, 2020, 1))
facts_west_2011_2020 <- filter(facts_west, YEAR %in% seq(2011, 2020, 1))

facts_2011_2020_U <- facts_union_yr(facts_2011_2020)
facts_west_2011_2020_U <- facts_union_yr(facts_west_2011_2020)

facts_2011_2020_U$Footprint_Acres <- as.numeric(st_area(facts_2011_2020_U)/4046.86)
facts_west_2011_2020_U$Footprint_Acres <- as.numeric(st_area(facts_west_2011_2020_U)/4046.86)

Footprint_Acres_2011_2020 <- sum(facts_2011_2020_U$Footprint_Acres)
Footprint_Acres_West_2011_2020 <- sum(facts_west_2011_2020_U$Footprint_Acres)

USFS_Budget <- read_excel(here("data", "raw", "USFS", "USFS_Budget_Hoover.xlsx"), skip = 1)

Cost_2011_2017 <- as.numeric(USFS_Budget[USFS_Budget$`Account/Activity` == "Hazardous Fuelse", c("FY2011", "FY2012", "FY2013", "FY2014",
                                                                                                 "FY2015", "FY2016", "FY2017")])

Cost_2018_2020 <- as.numeric(USFS_Budget[USFS_Budget$`Account/Activity` == "NFS Hazardous Fuelse", c("FY2018", "FY2019", "FY2020")])

TotalCost <- (sum(Cost_2011_2017) + sum(Cost_2018_2020))*1000000 # acount for millions

CPI_2020 <- as.numeric(filter(avg_cpi, START_YEAR == 2020)$CPI_FACTOR) # get 2020 to 2023 conversion factor 

TotalCost_Adj <- TotalCost / CPI_2020 # adjust to 2023 dollars 

# Average Cost per acre

TotalCost_Adj/Footprint_Acres_2011_2020 # will use in Cost-Benefit Analysis in script "SurvivalPlots.R" - $400.9842 per footprint acre


#### Back to code  ####

treatment_long <- yearly_acres_treated_USFS %>%
  rename(YEAR_MTBS = YEAR) %>%
  pivot_longer(cols = c(YEARLY_COST, ACRES_TREATED),
               names_to = "Metric",
               values_to = "Total") %>%
  mutate(Region = "Western USFS",
         Metric = recode(Metric,
                         "ACRES_TREATED" = "Acres of Fuel Treatment",
                         "YEARLY_COST" = "Fuel Treatment Costs"))


combined_data <- bind_rows(mtbs_long2, treatment_long)

# sets the order in which the plots are shown
combined_data$Metric <- factor(combined_data$Metric, levels = c(
  "Acres Burned", "Suppression Costs", "CO₂ Emissions", "PM₂.₅ Emissions",
  "Acres of Fuel Treatment","Fuel Treatment Costs", "Structures Destroyed"
))

# Get list of Metrics to plot
metric_list <- levels(combined_data$Metric)

# Create individual plots with lines only and shared theme
plot_list <- lapply(metric_list, function(metric_name) {
  ggplot(subset(combined_data, Metric == metric_name),
         aes(x = YEAR_MTBS, y = Total, color = Region, group = Region)) +
    geom_line(size = 1.2) +  # Thicker line for visibility
    scale_y_continuous(labels = scales::comma) +
    scale_color_manual(values = c("Western All" = "#E69F00", "Western USFS" = "#56B4E9")) +  # Subtle color scheme
    labs(
      title = metric_name,
      x = "Year",
      y = NULL,
      color = NULL
    ) +
    theme_minimal(base_size = 14) +  # Larger base size for readability
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Adjusted text size
      axis.text.y = element_text(size = 12),  # Make y-axis text more readable
      axis.title.x = element_text(size = 14, face = "bold"),  # Bold x-axis title
      axis.title.y = element_text(size = 14, face = "bold"),  # Bold y-axis title
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Bold and centered title
      legend.position = "none",  # Removed legend for all but last plot
      panel.grid.major = element_blank(),  # Removed major gridlines
      panel.grid.minor = element_blank(),  # Removed minor gridlines
      axis.line = element_line(size = 0.8, color = "black"),  # Add black axis lines
      plot.margin = margin(10, 10, 10, 10)  # Slight margin for clean edges
    )
})

# Add legend only to the last plot
plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
  theme(legend.position = "bottom")

# Arrange plots in a 2x4 grid
combined_plot <- patchwork::wrap_plots(plot_list, ncol = 4)

#### Map of Full Sample

## get full set of fires

MTBS_IDs <- unique(Grids_df$FIRE_ID)
mtbs_west <- filter(mtbs, WEST == 1)

mtbs_west <- mutate(mtbs_west, FT_FIRE_INT = ifelse(FIRE_ID %in% MTBS_IDs, "Yes", "No")) %>% 
  filter(YEAR_MTBS %in% seq(2017,2023, 1)) %>%
  st_transform(crs = 3857)

sum(mtbs_west$FT_FIRE_INT == "No")

as.data.frame(mtbs_west) %>%
  group_by(FT_FIRE_INT) %>%
  summarise(Total_Acres_Burned = sum(ACRES), 
            Total_Sup_Cost = sum(Total_Sup_Cost, na.rm = T))

burro_id <- filter(mtbs, FIRE_NAME == "BURRO" & YEAR_MTBS == 2017)$FIRE_ID

mtbs_fire_ig <- filter(mtbs_point, FIRE_ID == burro_id)  %>% st_transform(crs = 3857)

usa <- rnaturalearth::ne_states(country = "United States of America", returnclass = "sf") %>%
  st_transform(crs = 3857)

Western_States <- c("Washington", "Oregon", "California", "Idaho", "Utah", "Arizona",
                    "New Mexico", "Montana", "Colorado", "Wyoming", "Nevada")

usa_western <- filter(usa, woe_name %in% Western_States)

mtbs_west_1 <- sf::st_as_sf(mtbs_west) %>% 
  dplyr::filter(!is.na(FT_FIRE_INT)) %>%
  dplyr::mutate(
    FT_FIRE_INT = factor(as.character(FT_FIRE_INT),
                         levels = c("No", "Yes"))
  )

# High-quality map
tmap_mode("plot")  # Ensures static map for journal-quality output

# Expand the bounding box by adding a buffer (e.g., 0.05 units)
burro_buffer <- st_buffer(mtbs_fire_ig, dist = 50000) %>% st_transform(crs = 3857)

fire_FT_map <- tm_shape(usa_western) +
  tm_polygons(col = "grey90", border.col = "black", lwd = 1) +
  tm_shape(mtbs_west_1) +
  tm_polygons(
    col = "FT_FIRE_INT",                        # v3-compatible mapping
    palette = c("No" = "#D55E00", "Yes" = "#0072B2"),
    alpha = 0.7,
    border.col = "white", lwd = 0.3,
    
    # IMPORTANT: kill the auto legend for this layer
    legend.show = FALSE
  ) +
  # add legend manually (so no NA/Missing can appear)
  tm_add_legend(
    type = "fill",
    labels = c("No", "Yes"),
    col    = c("#D55E00", "#0072B2"),
    title  = "USFS fuel treatment interaction"
  ) +
  tm_layout(
    frame = FALSE,
    legend.outside = TRUE,
    fontfamily = "Helvetica",
    legend.position        = c("left", "bottom"),
    legend.stack = "horizontal",
    legend.title.size = 0.8,
    legend.text.size  = 0.7
  )

tmap_save(fire_FT_map, here("output", "maps", "SampleFires.png"),
          width = 12, height = 8, dpi = 300, outer.margins = 0)

img_path <- here("output", "maps", "SampleFires.png") 

trim_and_save <- function(input_path) { 
  img <- image_read(input_path) 
  img_trim <- image_trim(img) 
  image_write(img_trim, input_path) 
} 

trim_and_save(img_path)

map <- ggdraw() + draw_image(here("output", "maps", "SampleFires.png"))


combined_data_filt <- filter(combined_data, Metric != "PM₂.₅ Emissions" & Metric != "CO₂ Emissions" & 
                               Metric != "Structures Destroyed")

# sets the order in which the plots are shown
combined_data_filt$Metric <- factor(combined_data_filt$Metric, levels = c(
  "Acres Burned", "Suppression Costs",
  "Acres of Fuel Treatment","Fuel Treatment Costs"
))


# Get list of Metrics to plot
metric_list <- levels(combined_data_filt$Metric)

acres_FT <- combined_data_filt[combined_data_filt$Metric == "Fuel Treatment Costs" | combined_data_filt$Metric == "Acres of Fuel Treatment",]

acres_FT <- mutate(acres_FT, Region = "Western All")

combined_data_filt <- rbind(combined_data_filt, acres_FT)

combined_data_filt <- mutate(combined_data_filt, Region = ifelse(Region == "Western All", "Western U.S.", "USFS Western U.S."))

combined_data_filt$Region <- factor(combined_data_filt$Region,
                                    levels = c("Western U.S.", "USFS Western U.S."))

combined_data_filt$Total <- combined_data_filt$Total/1000000

# Create individual plots with lines only and shared theme
plot_list <- lapply(metric_list, function(metric_name) {
  
  y_lab <- dplyr::case_when(
    metric_name %in% c("Acres Burned", "Acres of Fuel Treatment") ~ "Million acres",
    metric_name %in% c("Suppression Costs", "Fuel Treatment Costs") ~ "Million US$ (2023)",
    TRUE ~ "Millions"
  )
  
  ggplot(subset(combined_data_filt, Metric == metric_name),
         aes(x = YEAR_MTBS, y = Total, color = Region, group = Region)) +
    geom_line(size = 1.2) +
    scale_y_continuous(labels = scales::comma) +
    scale_color_manual(
      values = c("Western U.S." = "#E69F00", "USFS Western U.S." = "#56B4E9"),
      drop = FALSE
    ) +
    labs(
      title = metric_name,
      x = "Year",
      y = y_lab,
      color = NULL
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y  = element_text(size = 7),
      axis.title.x = element_text(size = 8, face = "bold"),
      axis.title.y = element_text(size = 8, face = "bold"),
      plot.title   = element_text(size = 8.5, face = "bold", hjust = 0.5),
      legend.position   = "none",
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.line         = element_line(size = 0.8, color = "black"),
      plot.margin       = margin(10, 10, 10, 10)
    )
})

# Add legend only to the last plot
plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
  theme(legend.position = "bottom")

combined_plot <- wrap_plots(plot_list, ncol = 2) +
  plot_layout(guides = "collect") & 
  theme(
    legend.title = element_text(
      family = "Helvetica",
      face   = "bold",
      size   = 7
    ),
    legend.text = element_text(
      family = "Helvetica",
      size   = 7
    ),
    legend.position = "bottom"
  )


final_plot <- (map | combined_plot) + 
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.title = element_text(margin = margin(b = 1)),
        plot.tag = element_text(
          family = "Helvetica",
          face   = "bold",
          size   = 10
        ),
        panel.spacing = unit(0.8, "lines"),
        legend.position = "bottom",
        legend.title = element_text(family = "Helvetica", face = "bold", size = 7),
        legend.text  = element_text(family = "Helvetica", size = 7)
  )

# ggsave(here("output", "figures", "Figure1.pdf"), final_plot, width = 14, height = 10, units = "in")

ggsave(here("output", "figures", "Figure1.pdf"), final_plot, 
       width  = 7.24,          # Science 3-column width
       height = 5.5,    # ~4.5 in
       units  = "in",
       dpi    = 300)
