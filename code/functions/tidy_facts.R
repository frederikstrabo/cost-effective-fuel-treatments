##########################################    Tidy FACTS    ##########################################

################# Purpose: Cleans the full FACTS dataset, renames variables to match metadata,
#################      removes duplicates & unreliable observations, creates new variables, and adjusts for inflation

tidy_facts <- function (fs_fuel_trt, inf_yr, crs){
        
        fs_fuel_trt <- fs_fuel_trt[st_is_valid(fs_fuel_trt), ] # Get rid of invalid polygons
        
        fs_fuel_trt <- rename(fs_fuel_trt, ACTIVITY_CODE = ACTIVITY_C, LOCAL_QUALIFIER = LOCAL_QUAL, ASU_NBR_UNITS = ASU_NBR_UN, ADMIN_REGION_CODE = ADMIN_REGI, 
                          ADMIN_FOREST_CODE = ADMIN_FORE, ADMIN_DISTRICT_CODE = ADMIN_DIST, OWNERSHIP_CODE = OWNERSHIP_, PROC_REGION_CODE = PROC_REGIO,
                          PROC_FOREST_CODE = PROC_FORES, LAND_SUITABILITY_CLASS_CODE = LAND_SUITA, PRODUCTIVITY_CLASS_CODE = PRODUCTIVI, MGT_AREA_CODE = MGT_AREA_C,
                          MGT_PRESCRIPTION_CODE = MGT_PRESCR, NBR_UNITS_PLANNED = NBR_UNITS_, NBR_UNITS_ACCOMPLISHED = NBR_UNITS1, DATE_PLANNED = DATE_PLANN,
                          DATE_AWARDED = DATE_AWARD, DATE_COMPLETED = DATE_COMPL, FISCAL_YEAR_PLANNED = FISCAL_YEA, FISCAL_YEAR_COMPLETED = FISCAL_Y_1, 
                          FY_PLANNED_OR_ACCOMPLISHED = FY_PLANNED, METHOD_CODE = METHOD_COD, EQUIPMENT_CODE = EQUIPMENT_, COST_PER_UOM = COST_PER_U,
                          NEPA_PROJECT_ID = NEPA_PROJE, NEPA_DOC_NAME = NEPA_DOC_N, IMPLEMENTATION_PROJECT = IMPLEMENTA, IMPLEMENTATION_PROJECT_NBR = IMPLEMEN_1,
                          IMPLEMENTATION_PROJECT_TYPE = IMPLEMEN_2, ACCOMPLISHED_UNDER_HFRA = ACCOMPLISH, ACCOMPLISHED_UNDER_HFI = ACCOMPLI_1, 
                          ACTIVITY_CN = ACTIVITY_1, ACTIVITY_UNIT_CN = ACTIVITY_U, FEATURE_TYPE = FEATURE_TY, TREATMENT_TYPE = TREATMENT_, 
                          ACTIVITY_UNIT_NAME = ACTIVITY_2, ACTIVITY_SUB_UNIT_NAME = ACTIVITY_S, WORKFORCE_CODE = WORKFORCE_, NEPA_PROJECT_CN = NEPA_PRO_1,
                          TREATMENT_NAME = TREATMENT1, STAGE_VALUE = STAGE_VALU, DATA_SOURCE = DATA_SOURC, DATA_SOURCE_VALUE = DATA_SOU_1, 
                          FS_UNIT_NAME = FS_UNIT_NA, ETL_MODIFIED = ETL_MODIFI, EDW_INSERT_DATE = EDW_INSERT, ETL_MODIFIED_DATE_HAZ = ETL_MODI_1, 
                          ACT_CREATED_DATE = ACT_CREATE, ACT_MODIFIED_DATE = ACT_MODIFI)
        
        # Only look at Mechanical & RX Burns
        
        fs_fuel_trt <- rename(fs_fuel_trt, TREATMENT_CAT = CAT_NM) %>%
          filter(TREATMENT_CAT != "Other" & TREATMENT_CAT != "Wildfire Non-Treatment")
        
        fs_fuel_trt <- mutate(fs_fuel_trt, COMPLETE = ifelse(is.na(DATE_COMPLETED) == F, 1, 0))
        
        
        fs_fuel_trt_df <- as.data.frame(fs_fuel_trt)
        
        fs_fuel_trt_df$YEAR <-as.numeric(format(as.Date(fs_fuel_trt_df$DATE_COMPLETED) , "%Y")) 
        fs_fuel_trt_df$YEAR_PLANNED <-as.numeric(format(as.Date(fs_fuel_trt_df$DATE_PLANNED) , "%Y")) 
        fs_fuel_trt_df$MONTH <-as.numeric(format(as.Date(fs_fuel_trt_df$DATE_COMPLETED) , "%m")) 
        fs_fuel_trt_df$MONTH_PLANNED <-as.numeric(format(as.Date(fs_fuel_trt_df$DATE_PLANNED) , "%m")) 
        
        fs_fuel_trt_df_comp <- filter(fs_fuel_trt_df, COMPLETE == 1)
        fs_fuel_trt_df_incomp <- filter(fs_fuel_trt_df, COMPLETE == 0)
        
        # Inflation adjustment
        # library(quantmod)
        # getSymbols("CPIAUCSL", src='FRED') 
        cpi_data <- read_csv(here("data", "raw", "FRED", "CPIAUCSL.csv")) %>%
          mutate(YEAR = year(observation_date))
        
        # Average CPI by year
        avg_cpi <- cpi_data %>%
          group_by(YEAR) %>%
          summarise(AVG_CPI = mean(CPIAUCSL, na.rm = TRUE))
        
        # Reference CPI (e.g., 2020)
        ref_cpi <- avg_cpi %>%
          filter(YEAR == inf_yr) %>%
          pull(as.numeric(AVG_CPI))
        
        # Inflation adjustment factor
        avg_cpi <- avg_cpi %>%
          mutate(CPI_FACTOR = AVG_CPI / ref_cpi) %>%
          dplyr::select(YEAR, CPI_FACTOR)
        
        fs_fuel_trt_df_comp <- filter(fs_fuel_trt_df_comp, YEAR <= 2023) %>%
          left_join(avg_cpi, by = "YEAR") %>%
          mutate(ADJUSTED_COST_ACRE = COST_PER_UOM/CPI_FACTOR) %>%
          mutate(ADJUSTED_TOTAL_COST = NBR_UNITS_ACCOMPLISHED*ADJUSTED_COST_ACRE) %>%
          subset(select = -CPI_FACTOR)
        
        fs_fuel_trt_df_incomp$ADJUSTED_COST_ACRE <- NA
        fs_fuel_trt_df_incomp$ADJUSTED_TOTAL_COST <- NA
  
        fs_fuel_trt_df <- rbind(fs_fuel_trt_df_comp, fs_fuel_trt_df_incomp)
        
        fs_fuel_trt <- st_as_sf(fs_fuel_trt_df)
        fs_fuel_trt$Project_Acres <- as.numeric(st_area(fs_fuel_trt)/4046.86)

        
        fs_fuel_trt <- st_transform(fs_fuel_trt, crs = crs)
        fs_fuel_trt <- st_make_valid(fs_fuel_trt)
        
        return(fs_fuel_trt)
  
}

