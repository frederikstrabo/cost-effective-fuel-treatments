##########################################  02_cumulative_effects_matching.R  ##########################################

################# Purpose: Re-run the survival and BCR analysis using the matched control specification.

################# Outputs: Table S10 saved as "TableS10.tex" saved in "output/tables". 

################# Estimated run time: ~22 min.

rm(list=ls())

# if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr,sf, tmap, magrittr, rnaturalearth, rnaturalearthdata, ggplot2, maps, lwgeom, raster, stars, haven, stargazer, quantmod, lubridate, tidyr, ggpubr, 
               exactextractr, tictoc, terra, gtools, here, fixest, modelsummary, readr, rdrobust, prism, parallel,tmaptools, 
               OpenStreetMap, maptiles, gifski, geosphere, did, cowplot, didimputation, boot, survival, survminer, patchwork, gridExtra, Matching, cobalt, did, patchwork)

# Set Path
here::i_am("code/07_robustness/02_cumulative_effects_matching.R")

tic()

# Load Functions 

source(here("code", "functions","tidy_facts.R"))
source(here("code", "functions","calculate_size_survival.R"))
st_erase = function(x, y) st_difference(x, st_make_valid(st_union(st_combine(y)))) # Taken from here https://r-spatial.github.io/sf/reference/geos_binary_ops.html

#### Using code that does the matching from "07_robustness/04_conditional_effects_robustness.R"

###### Build Matched Subsample

Grids_df <- read_csv(here("data", "intermediate", "SpatialDiD_Grids_L24_K05.csv"))
Grids_df <- mutate(Grids_df, Treat_Post = ifelse(Treated_Dir == 1 & time_to_treat >= 0, 1, 0))
Grids_df <- mutate(Grids_df,  dist_treated = ifelse(Treated_Dir == 0, 1000, L_TAU)) # For Sun & Abraham Estimator

Grids_df_treated <- filter(Grids_df, time_to_treat > 0) %>% filter(Treated == 1 | Treated_LAG == 1) # For any observations post-treatment make sure they are still inside of a treatment or directly after one
Grids_df_yet_treated <- filter(Grids_df, time_to_treat <= 0)

Grids_df <- rbind(Grids_df_treated, Grids_df_yet_treated)

# subset to observations that are "yet-to-be extinguished" in order to estimate the conditional hazard model

Grids_df_filt <- filter(Grids_df, BURN_LAG == 1)

# limit to the 2.5km window - also remove treated observations that occur in distance bin 1 (i.e. Burn_Right_Away == 1).

Grids_df_filt_1 <- filter(Grids_df_filt, time_to_treat %in% seq(-5, 4,1) & Burn_Right_Away == 0)

# remove observations with no Delta MTT prediction for matching purposes (it errors if so)

Grids_df_filt_1 <- filter(Grids_df_filt_1, is.na(DELTA_MTT) == F)

set.seed(2020)

## X is the matrix of controls that we'd like to create balance on 

X  <- Grids_df_filt_1 %>% dplyr::select(distance_bin, USFS_NF, WindSpeed, ERC,
                                        DELTA_MTT, LOG_FIRE_INTENSITY, LOG_DIST_LAT, TRI)

Tr  <- Grids_df_filt_1$Treated_Dir == 1

Y  <- Grids_df_filt_1$BURN

# Matrix of controls to match on

BalanceMat <- cbind(Grids_df_filt_1$distance_bin, Grids_df_filt_1$USFS_NF, Grids_df_filt_1$WindSpeed, Grids_df_filt_1$ERC,
                    Grids_df_filt_1$DELTA_MTT, Grids_df_filt_1$LOG_FIRE_INTENSITY, Grids_df_filt_1$LOG_DIST_LAT, Grids_df_filt_1$TRI)

## Match observations using the GenMatch function from the "Matching" package

## Will match plots exactly based on their distance bin and whether or not it occurs in a National Forest or not. Plots will then inexactly matched to find the optimal 
##  covariate balance across the most important determinants of fire spread: wind speed, ERC, Delta T, log fire intensity (from MTT), 
##    log distance to nearest LAT drop, and topographic ruggedness index.

gen1 <- GenMatch(Tr = Tr, X = X,
                 BalanceMatrix = BalanceMat,
                 estimand = "ATT",
                 M = 1,
                 replace=FALSE,
                 exact = c(T, T, F, F, F, F, F, F),
                 pop.size=10,
                 max.generations=10,
                 wait.generations=1,
                 caliper = 0.15  # Define a caliper
)

mgens <- Match(Y=Y, Tr= Tr, X = X,
               estimand = "ATT",
               M = 1,
               replace=FALSE,
               exact = c(T, T, F, F, F, F, F, F),
               Weight.matrix = gen1)

balance.table_tr <- bal.tab(mgens, Treated_Dir ~
                              USFS_NF + WindSpeed + ERC + DELTA_MTT + LOG_FIRE_INTENSITY + LOG_DIST_LAT + TRI, data = Grids_df_filt_1,
                            un = TRUE)

## subset to matched sample

data.matched <- bind_rows(Grids_df_filt_1 %>% slice(mgens$index.treated), # subset data
                          Grids_df_filt_1 %>% slice(mgens$index.control))

###### Now that we have our matched sample we will use code from "06_analysis/03_cumulative_effects.R" to do survival and BCR analysis

#### Predict out of sample

controls <- c("Slope", "Elev", "TRI", "Distance_FS_Road", "USFS_NF", "Wilderness", "MFRI", 
              "Distance_US_Highway", "Distance_WUI", "ERC", "WindSpeed", "FM1000", "PREV_BURN_10Y",
              "LAT", "LAT_LINE_INT", "LOG_DIST_LAT", "WIND_DIFF", "DELTA_MTT", "NA_DELTA_MTT",
              "LOG_FIRE_INTENSITY", "ROAD", "WUI")

# Dynamically construct the formula
formula_burn_pred <- as.formula(paste(
  "BURN ~ ",
  paste(controls, collapse = " + "),           # Add all control variables
  "| distance_bin + direction_fire_FE + FuelType_2001 + Aspect_Class"         # Add fixed effects
))


# Dynamically construct the formula
formula_burn_sev <- as.formula(paste(
  "BURN_SEV0 ~ ",
  paste(controls, collapse = " + "),           # Add all control variables
  "| distance_bin + direction_fire_FE + FuelType_2001 + Aspect_Class"         # Add fixed effects
))

# Dynamically construct the formula
formula_med_high_burn_sev <- as.formula(paste(
  "BURN_SEV_MED_HIGH_PCT ~ ",
  paste(controls, collapse = " + "),           # Add all control variables
  "| direction_fire_FE + distance_bin + FuelType_2001 + Aspect_Class"         # Add fixed effects
))

formula_high_burn_sev <- as.formula(paste(
  "BURN_SEV_HIGH_PCT ~ ",
  paste(controls, collapse = " + "),           # Add all control variables
  "| direction_fire_FE + distance_bin + FuelType_2001 + Aspect_Class"         # Add fixed effects
))

## Create untreated samples - post-treatment observations are considered treated if they are inside of a treatment (Treated) or directly after one (Treated_LAG) 

Grids_df_untreated <- filter(data.matched, Treated == 0 & Treated_LAG == 0)
Grids_df_untreated_BURN <- filter(Grids_df_untreated, BURN == 1)

## Run OLS models

OLS_Pred_Model <- feols(formula_burn_pred,
                        cluster = ~FIRE_ID,
                        weights = ~Grid_Acres,
                        fixef.rm  = "none",
                        data = Grids_df_untreated)

OLS_Pred_Model_Burn <- feols(formula_burn_sev,
                             cluster = ~FIRE_ID,
                             weights = ~Grid_Acres,
                             fixef.rm  = "none",
                             data = Grids_df_untreated_BURN)  

OLS_Pred_Model_Med_High_Burn <- feols(formula_med_high_burn_sev,
                                      cluster = ~FIRE_ID,
                                      weights = ~Grid_Acres,
                                      fixef.rm  = "none",
                                      data = Grids_df_untreated_BURN)

OLS_Pred_Model_High_Burn <- feols(formula_high_burn_sev,
                                  cluster = ~FIRE_ID,
                                  weights = ~Grid_Acres,
                                  fixef.rm  = "none",
                                  data = Grids_df_untreated_BURN)

#### Define distance  window that we allow cumulative effects to accumulate in B-C Analysis

distance_window <- seq(-1,9,1)

treated_data <- Grids_df %>% filter(Treated_Dir == 1 & Burn_Right_Away == 0) %>% filter(time_to_treat %in% distance_window)

## Predict conditional probabilities & burn severities

treated_data$pred_without_treatment <- predict(OLS_Pred_Model, newdata = treated_data, type = "response")
treated_data$pred_without_treatment_BS <- predict(OLS_Pred_Model_Burn, newdata = treated_data, type = "response")
treated_data$pred_without_treatment_MED_HIGH_BS <- predict(OLS_Pred_Model_Med_High_Burn, newdata = treated_data, type = "response")
treated_data$pred_without_treatment_HIGH_BS <- predict(OLS_Pred_Model_High_Burn, newdata = treated_data, type = "response")


sum(is.na(treated_data$pred_without_treatment))
sum(is.na(treated_data$pred_without_treatment_BS))
sum(is.na(treated_data$pred_without_treatment_MED_HIGH_BS))
sum(is.na(treated_data$pred_without_treatment_HIGH_BS))

## Remove observations for which we have no prediction

treated_data <- filter(treated_data, is.na(pred_without_treatment) == F & is.na(pred_without_treatment_BS) == F)

sum(is.na(treated_data$pred_without_treatment))

ex <- dplyr::select(treated_data, direction_fire_FE, BURN, time_to_treat, dist_treated, distance_bin, pred_without_treatment) %>% arrange(direction_fire_FE)

treated_data[is.na(treated_data$pred_without_treatment), "pred_without_treatment"] <- 0

## Option 1 - for predictions > 1 set to 1, < 0, set to 0

treated_data <- mutate(treated_data, pred_without_treatment = case_when(pred_without_treatment > 1 ~ 1,
                                                                        pred_without_treatment < 0 ~ 0,
                                                                        pred_without_treatment >= 0 & pred_without_treatment <= 1 ~ pred_without_treatment),
                       pred_without_treatment_MED_HIGH_BS = case_when(pred_without_treatment_MED_HIGH_BS > 1 ~ 1,
                                                                      pred_without_treatment_MED_HIGH_BS < 0 ~ 0,
                                                                      pred_without_treatment_MED_HIGH_BS >= 0 & pred_without_treatment_MED_HIGH_BS <= 1 ~ pred_without_treatment_MED_HIGH_BS),
                       pred_without_treatment_HIGH_BS = case_when(pred_without_treatment_HIGH_BS > 1 ~ 1,
                                                                  pred_without_treatment_HIGH_BS < 0 ~ 0,
                                                                  pred_without_treatment_HIGH_BS >= 0 & pred_without_treatment_HIGH_BS <= 1 ~ pred_without_treatment_HIGH_BS)
)


treated_data[is.na(treated_data$pred_without_treatment_BS), "pred_without_treatment_BS"] <- 0
treated_data[is.na(treated_data$pred_without_treatment_MED_HIGH_BS), "pred_without_treatment_MED_HIGH_BS"] <- 0
treated_data[is.na(treated_data$pred_without_treatment_HIGH_BS), "pred_without_treatment_HIGH_BS"] <- 0

# treated_data <- mutate(treated_data, pred_without_treatment_HIGH_BS = case_when(pred_without_treatment_HIGH_BS > 1 ~ 1,
#                                                                                 pred_without_treatment_HIGH_BS < 0 ~ 0,
#                                                                                 pred_without_treatment_HIGH_BS >= 0 & pred_without_treatment_HIGH_BS <= 1 ~ pred_without_treatment_HIGH_BS))

## "With treatment" counterfactual is observed outcomes

treated_data <- mutate(treated_data, 
                       pred_with_treatment = BURN,
                       pred_with_treatment_BS = BURN_SEV0,
                       pred_with_treatment_MED_HIGH_BS = BURN_SEV_MED_HIGH_PCT,
                       pred_with_treatment_HIGH_BS = BURN_SEV_HIGH_PCT) %>%
  arrange(direction_fire_FE, time_to_treat)


# Create 1-9 lags for the conditional probability of fire spread
treated_data <- treated_data %>%
  arrange(direction_fire_FE, time_to_treat) %>% # Ensure data is ordered by ID and time
  group_by(direction_fire_FE) %>%              # Group by ID
  mutate(p1 = lag(pred_with_treatment, 1),   
         p2 = lag(pred_with_treatment, 2),
         p3 = lag(pred_with_treatment, 3),
         p4 = lag(pred_with_treatment, 4),
         p5 = lag(pred_with_treatment, 5),
         p6 = lag(pred_with_treatment, 6),
         p7 = lag(pred_with_treatment, 7),
         p8 = lag(pred_with_treatment, 8),
         p9 = lag(pred_with_treatment, 9),
         p100 = lag(pred_without_treatment, 1),
         p200 = lag(pred_without_treatment, 2),
         p300 = lag(pred_without_treatment, 3),
         p400 = lag(pred_without_treatment, 4),
         p500 = lag(pred_without_treatment, 5),
         p600 = lag(pred_without_treatment, 6),
         p700 = lag(pred_without_treatment, 7),
         p800 = lag(pred_without_treatment, 8),
         p900 = lag(pred_without_treatment, 9)) %>% # Create lagged columns
  ungroup()                                      # Remove grouping

## Create survival probabilities - the production of conditional probabilities

treated_data <- mutate(treated_data, 
                       survival_with_treatment = case_when(time_to_treat == -1 ~ 1,
                                                           time_to_treat == 0 ~ pred_with_treatment,
                                                           time_to_treat == 1 ~ pred_with_treatment*p1, 
                                                           time_to_treat == 2 ~ pred_with_treatment*p1*p2,
                                                           time_to_treat == 3 ~ pred_with_treatment*p1*p2*p3,
                                                           time_to_treat == 4 ~ pred_with_treatment*p1*p2*p3*p4,
                                                           time_to_treat == 5 ~ pred_with_treatment*p1*p2*p3*p4*p5,
                                                           time_to_treat == 6 ~ pred_with_treatment*p1*p2*p3*p4*p5*p6,
                                                           time_to_treat == 7 ~ pred_with_treatment*p1*p2*p3*p4*p5*p6*p7,
                                                           time_to_treat == 8 ~ pred_with_treatment*p1*p2*p3*p4*p5*p6*p7*p8,
                                                           time_to_treat == 9 ~ pred_with_treatment*p1*p2*p3*p4*p5*p6*p7*p8*p9),
                       
                       survival_without_treatment = case_when(time_to_treat == -1 ~ 1,
                                                              time_to_treat == 0 ~ pred_without_treatment,
                                                              time_to_treat == 1 ~ pred_without_treatment*p100, 
                                                              time_to_treat == 2 ~ pred_without_treatment*p100*p200,
                                                              time_to_treat == 3 ~ pred_without_treatment*p100*p200*p300,
                                                              time_to_treat == 4 ~ pred_without_treatment*p100*p200*p300*p400,
                                                              time_to_treat == 5 ~ pred_without_treatment*p100*p200*p300*p400*p500,
                                                              time_to_treat == 6 ~ pred_without_treatment*p100*p200*p300*p400*p500*p600,
                                                              time_to_treat == 7 ~ pred_without_treatment*p100*p200*p300*p400*p500*p600*p700,
                                                              time_to_treat == 8 ~ pred_without_treatment*p100*p200*p300*p400*p500*p600*p700*p800,
                                                              time_to_treat == 9 ~ pred_without_treatment*p100*p200*p300*p400*p500*p600*p700*p800*p900)
)

## Create burn severity predicted probabilities product of the predicted survival probability * predicted conditional burn severity

treated_data <- mutate(treated_data, 
                       BS_with_treatment = survival_with_treatment*pred_with_treatment_BS,
                       BS_without_treatment = ifelse(time_to_treat == -1, survival_with_treatment*pred_with_treatment_BS,
                                                     survival_without_treatment*pred_without_treatment_BS),
                       MED_HIGH_BS_with_treatment =  survival_with_treatment*pred_with_treatment_MED_HIGH_BS,
                       MED_HIGH_BS_without_treatment = ifelse(time_to_treat == -1, survival_with_treatment*pred_with_treatment_MED_HIGH_BS,
                                                              survival_without_treatment*pred_without_treatment_MED_HIGH_BS),
                       HIGH_BS_with_treatment =  survival_with_treatment*pred_with_treatment_HIGH_BS,
                       HIGH_BS_without_treatment = ifelse(time_to_treat == -1, survival_with_treatment*pred_with_treatment_HIGH_BS,
                                                          survival_without_treatment*pred_without_treatment_HIGH_BS)
)

## Note how many directions there are

N_obs_treat_time <- treated_data %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n())

N_Directions <- as.numeric(N_obs_treat_time[1,2])

## Calculate average outcomes in the treated and untreated counterfactuals for each event time interval

survivor_plot <- treated_data %>%
  group_by(time_to_treat) %>%
  summarise(Mean_Survival_Treat = (sum(survival_with_treatment, na.rm = T)/N_Directions),
            Mean_Survival_No_Treat = (sum(survival_without_treatment, na.rm = T)/N_Directions),
            Mean_BS_Treat = sum(BS_with_treatment, na.rm = T)/N_Directions,
            Mean_BS_No_Treat = sum(BS_without_treatment, na.rm = T)/N_Directions,
            Pct_Burned_MS_HS_Treat = sum(MED_HIGH_BS_with_treatment, na.rm = T)/N_Directions,
            Pct_Burned_MS_HS_No_Treat = sum(MED_HIGH_BS_without_treatment, na.rm = T)/N_Directions,
            Pct_Burned_HS_Treat = sum(HIGH_BS_with_treatment, na.rm = T)/N_Directions,
            Pct_Burned_HS_No_Treat = sum(HIGH_BS_without_treatment, na.rm = T)/N_Directions,
            Diff_Survival = Mean_Survival_Treat - Mean_Survival_No_Treat,
            Diff_BS = Mean_BS_Treat - Mean_BS_No_Treat,
            Diff_MED_HIGH_BS = Mean_BS_Treat - Mean_BS_No_Treat,
            Diff_HIGH_BS = Pct_Burned_HS_Treat - Pct_Burned_HS_No_Treat,
            Acres_Burned_MS_HS_Treat = sum(MED_HIGH_BS_with_treatment*Grid_Acres),
            Acres_Burned_MS_HS_No_Treat = sum(MED_HIGH_BS_without_treatment*Grid_Acres),
            Acres_Burned_HS_Treat = sum(HIGH_BS_with_treatment*Grid_Acres),
            Acres_Burned_HS_No_Treat = sum(HIGH_BS_without_treatment*Grid_Acres)) %>%
  filter(time_to_treat %in% seq(-1,4,1))

survivor_plot <- mutate(survivor_plot, Pct_BS_Change = (Mean_BS_No_Treat - Mean_BS_Treat)/Mean_BS_No_Treat)

survivor_plot$distance <- seq(0, 2.5, 0.5)

interval_labels <- c(
  "0", "0.5", "1",
  "1.5", "2", "2.5"
)


# Plotting the time series for predicted survival rates
survival_plot <- ggplot(survivor_plot, aes(x = time_to_treat)) +
  # Dashed vertical line at -0.5
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_line(aes(y = Mean_Survival_Treat, color = "With Treatment (Observed)"), size = 1.2) +
  geom_line(aes(y = Mean_Survival_No_Treat, color = "Without Treatment (Counterfactual)"), size = 1.2) +
  labs(
    title = "Burn Probability by Distance From Treatment",
    x = "Distance from treatment interaction (0.5km)",
    y = "Mean Predicted Burn Probability",
    color = "Treatment Status"
  ) +
  scale_color_manual(values = c("With Treatment (Observed)" = "#0072B2", "Without Treatment (Counterfactual)" = "#D55E00")) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),  # Remove minor gridlines for a cleaner look
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add axes
    axis.line = element_blank(),  # Avoid duplicate axis lines
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

survival_plot


##################    Benefit-Cost Analysis   ##################

#### Counterfactual Acres Burned - With and Without Treatment

## Using survival probabilities multiply it by the acres (or assets at risk) in a plot to estimate acres burned (or other savings) under the no treatment counterfactual

treated_data <- mutate(treated_data, Acres_Burned_With_Treat = survival_with_treatment*Grid_Acres,
                       Acres_Burned_No_Treat = survival_without_treatment*Grid_Acres,
                       Acres_Burned_Realized = BURN*Grid_Acres,
                       Acres_Burned_MS_HS_With_Treat = survival_with_treatment*Grid_Acres*pred_with_treatment_MED_HIGH_BS, 
                       Acres_Burned_MS_HS_Without_Treat = survival_without_treatment*Grid_Acres*pred_without_treatment_MED_HIGH_BS,
                       Acres_Burned_MS_HS_Realized = BURN*Grid_Acres*BURN_SEV_MED_HIGH_PCT,
                       Acres_Burned_HS_With_Treat = survival_with_treatment*Grid_Acres*pred_with_treatment_HIGH_BS, 
                       Acres_Burned_HS_Without_Treat = survival_without_treatment*Grid_Acres*pred_without_treatment_HIGH_BS,
                       Acres_Burned_HS_Realized = BURN*Grid_Acres*BURN_SEV_HIGH_PCT
)

sum(treated_data$Acres_Burned_With_Treat)
sum(treated_data$Acres_Burned_No_Treat)
sum(treated_data$Acres_Burned_Realized)

acres_burned_savings <- sum(treated_data$Acres_Burned_No_Treat, na.rm = T) - sum(treated_data$Acres_Burned_With_Treat, na.rm = T) # Savings only considering 5 leads
acres_burned_savings_pct <- (sum(treated_data$Acres_Burned_No_Treat, na.rm = T) - sum(treated_data$Acres_Burned_With_Treat, na.rm = T))/sum(treated_data$Acres_Burned_Realized, na.rm = T)

acres_burned_savings
acres_burned_savings_pct


acres_burned_MS_HS_savings <- sum(treated_data$Acres_Burned_MS_HS_Without_Treat, na.rm = T) - sum(treated_data$Acres_Burned_MS_HS_With_Treat, na.rm = T) # Savings only considering 5 leads
acres_burned_MS_HS_savings_pct <- (sum(treated_data$Acres_Burned_MS_HS_Without_Treat, na.rm = T) - sum(treated_data$Acres_Burned_MS_HS_With_Treat, na.rm = T))/sum(treated_data$Acres_Burned_MS_HS_Realized)

acres_burned_MS_HS_savings
acres_burned_MS_HS_savings_pct

acres_burned_HS_savings <- sum(treated_data$Acres_Burned_HS_Without_Treat, na.rm = T) - sum(treated_data$Acres_Burned_HS_With_Treat, na.rm = T) # Savings only considering 5 leads
acres_burned_HS_savings_pct <- (sum(treated_data$Acres_Burned_HS_Without_Treat, na.rm = T) - sum(treated_data$Acres_Burned_HS_With_Treat, na.rm = T))/sum(treated_data$Acres_Burned_HS_Without_Treat, na.rm = T)

acres_burned_HS_savings
acres_burned_HS_savings_pct

treated_data <- mutate(treated_data, Structure_Burned_With_Treat = survival_with_treatment*Struc_Count,
                       Structure_Burned_No_Treat = survival_without_treatment*Struc_Count,
                       Structure_Burned_Realized = BURN*Struc_Count, 
                       Structure_Value_Lost_With_Treat = survival_with_treatment*TOT_STRUC_VAL,
                       Structure_Value_Lost_No_Treat = survival_without_treatment*TOT_STRUC_VAL,
                       Structure_Value_Lost_Realized = BURN*TOT_STRUC_VAL)

struc_savings <- sum(treated_data$Structure_Burned_No_Treat, na.rm = T) - sum(treated_data$Structure_Burned_With_Treat, na.rm = T) 
struc_savings # Savings only considering 5 leads


(sum(treated_data$Structure_Burned_No_Treat, na.rm = T) - sum(treated_data$Structure_Burned_With_Treat, na.rm = T))/sum(treated_data$Structure_Burned_Realized, na.rm = T)

#### Structure Value Savings

house_val_savings <- sum(treated_data$Structure_Value_Lost_No_Treat, na.rm = T) - sum(treated_data$Structure_Value_Lost_With_Treat, na.rm = T) 
house_val_savings

(sum(treated_data$Structure_Value_Lost_No_Treat, na.rm = T) - sum(treated_data$Structure_Value_Lost_With_Treat, na.rm = T))/sum(treated_data$Structure_Value_Lost_Realized, na.rm = T)

#### Smoke Emissions

treated_data <- mutate(treated_data, PM25_With_Treat = Acres_Burned_With_Treat*FIRE_PM25_PER_ACRE,
                       PM25_No_Treat = Acres_Burned_No_Treat*FIRE_PM25_PER_ACRE,
                       PM25_Realized = Acres_Burned_Realized*FIRE_PM25_PER_ACRE,
                       CO2_With_Treat = Acres_Burned_With_Treat*FIRE_CO2_PER_ACRE,
                       CO2_No_Treat = Acres_Burned_No_Treat*FIRE_CO2_PER_ACRE,
                       CO2_Realized = Acres_Burned_Realized*FIRE_CO2_PER_ACRE)


PM25_savings <- sum(treated_data$PM25_No_Treat, na.rm = T) - sum(treated_data$PM25_With_Treat, na.rm = T) # PM2.5 Savings
PM25_savings_pct <- (sum(treated_data$PM25_No_Treat, na.rm = T) - sum(treated_data$PM25_With_Treat, na.rm = T))/sum(treated_data$PM25_Realized, na.rm = T)

CO2_savings <- sum(treated_data$CO2_No_Treat, na.rm = T) - sum(treated_data$CO2_With_Treat, na.rm = T) # CO2 Savings
CO2_savings_pct <- (sum(treated_data$CO2_No_Treat, na.rm = T) - sum(treated_data$CO2_With_Treat, na.rm = T))/sum(treated_data$CO2_Realized, na.rm = T)

PM25_savings
CO2_savings



carbon_savings <- (sum(treated_data$CO2_No_Treat, na.rm = T) - sum(treated_data$CO2_With_Treat, na.rm = T))*185 # CO2 * Carbon price
carbon_savings


#### Health Costs

treated_data <- mutate(treated_data, 
                       Deaths_With_Treat = Acres_Burned_With_Treat*Deaths_PER_ACRE,
                       Deaths_No_Treat = Acres_Burned_No_Treat*Deaths_PER_ACRE,
                       Deaths_Realized = Acres_Burned_Realized*Deaths_PER_ACRE,
                       DeathCost_With_Treat = Acres_Burned_With_Treat*DeathCost_PER_ACRE,
                       DeathCost_No_Treat = Acres_Burned_No_Treat*DeathCost_PER_ACRE,
                       DeathCost_Realized = Acres_Burned_Realized*DeathCost_PER_ACRE,
                       EarningsLoss_With_Treat = Acres_Burned_With_Treat*EarningsLoss_PER_ACRE,
                       EarningsLoss_No_Treat = Acres_Burned_No_Treat*EarningsLoss_PER_ACRE,
                       EarningsLoss_Realized = Acres_Burned_Realized*EarningsLoss_PER_ACRE)


Saved_Deaths <- sum(treated_data$Deaths_No_Treat, na.rm = T) - sum(treated_data$Deaths_With_Treat, na.rm = T) # Saved Deaths
Saved_Deaths_pct <- (sum(treated_data$Deaths_No_Treat, na.rm = T) - sum(treated_data$Deaths_With_Treat, na.rm = T))/sum(treated_data$Deaths_Realized, na.rm = T)

Saved_Deaths

Death_savings <- (sum(treated_data$DeathCost_No_Treat, na.rm = T) - sum(treated_data$DeathCost_With_Treat, na.rm = T)) # Deaths*VSL
Death_savings_pct <- (sum(treated_data$DeathCost_No_Treat, na.rm = T) - sum(treated_data$DeathCost_With_Treat, na.rm = T))/sum(treated_data$DeathCost_Realized, na.rm = T)

Death_savings

Earnings_savings <- (sum(treated_data$EarningsLoss_No_Treat, na.rm = T) - sum(treated_data$EarningsLoss_With_Treat, na.rm = T)) # Earnings Estimate
Earnings_savings_pct <- (sum(treated_data$EarningsLoss_No_Treat, na.rm = T) - sum(treated_data$EarningsLoss_With_Treat, na.rm = T))/sum(treated_data$EarningsLoss_Realized, na.rm = T)

Earnings_savings

Health_Savings <- Death_savings + Earnings_savings


######## Calculate treatment costs for the treatments in our sample & Estimate the probability of treatment intersecting with a fire over its life time for all treatments in the Western U.S. ########

#### FACTS - Forest Service Fuel Treatments

facts <- st_read(here("data", "raw", "FACTS", "S_USA.Activity_HazFuelTrt_PL.shp"))
facts <- tidy_facts(facts, inf_yr = 2023, crs = 5070)
facts$ROW_ID <- 1:nrow(facts)

states <- c("WA", "OR", "CA", "ID", "NV", "AZ", "MT", "UT", "NM", "WY", "CO")

facts_west <- filter(facts, COMPLETE == 1) %>% filter(STATE_ABBR %in% states)

facts_west_filt <- filter(facts_west, YEAR >= 2004)

#### Load in MTBS

# Burn Area Boundaries

mtbs <- st_read(here("data", "raw", "MTBS", "MTBS_Burn_Area", "S_USA.MTBS_BURN_AREA_BOUNDARY.shp")) %>%
  st_transform(mtbs, crs = 5070) %>%
  st_make_valid() %>%
  dplyr::rename(YEAR_MTBS = YEAR, MONTH_MTBS = STARTMONTH)

mtbs <- filter(mtbs, FIRE_TYPE == "Unknown" | FIRE_TYPE == "Wildfire") # No Wildland Fire Use Fires

mtbs_filt <- filter(mtbs, YEAR_MTBS >= 2004)

######## Estimate the probability of fuel treatment interacting with a fire over its life time by using survival function approach - Kaplan-Meier Estimator

## For each project in FACTS take its union and note when its completed

facts_act_grouped <- facts_west_filt %>%
  group_by(ACTIVITY_UNIT_CN) %>%
  summarise(geometry = st_union(geometry), 
            YEAR_COMP = max(YEAR), 
            YEAR_COMP_MIN = min(YEAR))

# Calculate acres treated by its footprint
facts_act_grouped$Acres_Treated <- as.numeric(st_area(facts_act_grouped)/4046.86)

# intersect projects with MTBS
facts_mtbs_int <- st_intersection(facts_act_grouped, mtbs_filt)

# Only count projects that intersect with projects completed at least 10 years prior to fire
activities_int <- as.data.frame(facts_mtbs_int) %>%
  filter(YEAR_COMP <= YEAR_MTBS | YEAR_COMP_MIN <= YEAR_MTBS) %>%
  filter(YEAR_MTBS - YEAR_COMP <= 10 | YEAR_MTBS - YEAR_COMP_MIN <= 10) %>%
  group_by(ACTIVITY_UNIT_CN) %>%
  summarise(INTERSECTS = 1)

facts_act_grouped <- merge(facts_act_grouped, activities_int, by = "ACTIVITY_UNIT_CN", all.x = T)

facts_act_grouped[is.na(facts_act_grouped$INTERSECTS), "INTERSECTS"] <- 0

# Compute time-to-event and censoring indicator
facts_act_grouped$event_time <- pmin(2023 - facts_act_grouped$YEAR_COMP, 10)  # Max follow-up is 10 years

facts_act_grouped <- as.data.frame(facts_act_grouped)

#### Estimate effects by size category: small = 75-600 acres, medium = 600-2400 acres, large > 2,400 acres.

treated_data_filt <- filter(treated_data, YEAR %in% seq(2017,2023,1))
quantiles <- quantile(treated_data_filt$TREAT_SIZE, probs = c(1/3,2/3,1)) 

q1 <- 600
q2 <- 2400
q_05 <- quantile(treated_data_filt$TREAT_SIZE, probs = .05) 

facts_act_grouped_small <- filter(facts_act_grouped, Acres_Treated > q_05 & Acres_Treated < q1)
facts_act_grouped_medium <- filter(facts_act_grouped, Acres_Treated >= q1 & Acres_Treated <= q2)
facts_act_grouped_large <- filter(facts_act_grouped, Acres_Treated > q2)

# Fit Kaplan-Meier survival curve
km_fit_1 <- survfit(Surv(event_time, INTERSECTS) ~ 1, data = facts_act_grouped_small)
km_fit_2 <- survfit(Surv(event_time, INTERSECTS) ~ 1, data = facts_act_grouped_medium)
km_fit_3 <- survfit(Surv(event_time, INTERSECTS) ~ 1, data = facts_act_grouped_large)

# Estimate lambda (probability of intersection within 10 years)
lambda_1 <- 1 - summary(km_fit_1, times = 10)$surv
lambda_2 <- 1 - summary(km_fit_2, times = 10)$surv
lambda_3 <- 1 - summary(km_fit_3, times = 10)$surv

treated_data_filt <- mutate(treated_data_filt, lambda = case_when(TREAT_SIZE < q1 ~ lambda_1,
                                                                  TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ lambda_2,
                                                                  TREAT_SIZE > q2 ~ lambda_3,))

treated_data_filt <- mutate(treated_data_filt, Size_Bin = case_when(TREAT_SIZE < q1 ~ "Small",
                                                                    TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ "Medium",
                                                                    TREAT_SIZE > q2 ~ "Large",))

treat_ids <- unique(treated_data_filt$TREAT_ID)
treat_ids <- treat_ids[!is.na(treat_ids)]

treatment_costs <- treated_data_filt %>%
  group_by(TREAT_ID) %>%
  summarise(TREAT_COST = dplyr::first(TREAT_COST), 
            lambda = dplyr::first(lambda),
            Size_Bin = dplyr::first(Size_Bin),
            YEAR_COMPLETE = dplyr::first(YEAR) - dplyr::first(YEARS_SINCE_TREAT) ,
            TREAT_SIZE = dplyr::first(TREAT_SIZE)
  )


Treat_Cost <- sum(treatment_costs$TREAT_COST, na.rm = T)

## Identify observations with high cost/acre - replace such treatments with median cost-acre

mean_cost <- mean(facts_west$COST_PER_UOM, na.rm = TRUE)
sd_cost <- sd(facts_west$COST_PER_UOM, na.rm = TRUE)
threshold <- mean_cost + 10*sd_cost
median_cost <- median(facts_west$COST_PER_UOM, na.rm = TRUE)

outliers <- unique(filter(facts_west, COST_PER_UOM > threshold)$ACTIVITY_UNIT_CN)
treatment_costs <- mutate(treatment_costs, TREAT_COST_REP = ifelse(TREAT_ID %in% outliers, TREAT_SIZE*median_cost, TREAT_COST))

Treat_Cost <- sum(treatment_costs$TREAT_COST_REP, na.rm = T)
Avoided_Damages_Matching <- house_val_savings + carbon_savings + Health_Savings
Benefit_Cost_Ratio_Conditional = (Avoided_Damages_Matching)/Treat_Cost

###### First calculate average benefit, conditional on fire arrival, per size class 

## Aggregate savings in each size bin

Treat_Size_Benefits <- treated_data_filt %>%
  group_by(Size_Bin) %>%
  summarise(Acres_Burned_Savings = sum(Acres_Burned_No_Treat - Acres_Burned_With_Treat, na.rm = T),
            Struc_Savings = sum(Structure_Burned_No_Treat - Structure_Burned_With_Treat, na.rm = T),
            House_Savings = sum(Structure_Value_Lost_No_Treat - Structure_Value_Lost_With_Treat, na.rm = T),
            Carbon_Savings = (sum(CO2_No_Treat - CO2_With_Treat, na.rm = T))*185,
            Death_Savings = sum(DeathCost_No_Treat - DeathCost_With_Treat, na.rm = T),
            Earning_Savings = sum(EarningsLoss_No_Treat - EarningsLoss_With_Treat, na.rm = T),
            N_Treats = n_distinct(TREAT_ID),
            lambda = dplyr::first(lambda))

Treat_Size_Benefits <- mutate(Treat_Size_Benefits, Total_Benefit = House_Savings + Carbon_Savings + Death_Savings)
Treat_Size_Benefits <- mutate(Treat_Size_Benefits, Benefit_Per_Treat = Total_Benefit/N_Treats)

#### Benefit Cost Ratios for different sized treatments

mu_1 <- filter(Treat_Size_Benefits, Size_Bin == "Small")$Benefit_Per_Treat
mu_2 <- filter(Treat_Size_Benefits, Size_Bin == "Medium")$Benefit_Per_Treat
mu_3 <- filter(Treat_Size_Benefits, Size_Bin == "Large")$Benefit_Per_Treat

#### Calculate the Ex-ante B/C Ratio for USFS treatments that could have intersected with the wildfires in our 
####    sample during their effective lifetime (10 years) - i.e. those completed from 2007-2023. 

facts_act_grouped_06_23 <- facts_west %>%
  filter(YEAR >= 2007 & YEAR <= 2023) %>%
  group_by(ACTIVITY_UNIT_CN) %>%
  summarise(geometry = st_union(geometry), 
            TOT_COST = sum(ADJUSTED_TOTAL_COST, is.na = T),
            COST_PER_UOM = max(COST_PER_UOM),
            ACCOMPLISHED_ACRES = sum(NBR_UNITS_ACCOMPLISHED)
  )

facts_act_grouped_06_23$TREAT_SIZE <- as.numeric(st_area(facts_act_grouped_06_23)/4046.86)

facts_act_grouped_06_23 <- as.data.frame(facts_act_grouped_06_23)

## For each treatment note its estimated probability of fire interaction (lambda) and its benefit conditional on fire (mu) based on its size.

facts_act_grouped_06_23 <- mutate(facts_act_grouped_06_23, lambda = case_when(TREAT_SIZE < q1 ~ lambda_1,
                                                                              TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ lambda_2,
                                                                              TREAT_SIZE > q2 ~ lambda_3),
                                  mu = case_when(TREAT_SIZE < q1 ~ mu_1,
                                                 TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ mu_2,
                                                 TREAT_SIZE > q2 ~ mu_3))

facts_act_grouped_06_23 <- mutate(facts_act_grouped_06_23, Size_Bin = case_when(TREAT_SIZE < q1 ~ "Small",
                                                                                TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ "Medium",
                                                                                TREAT_SIZE > q2 ~ "Large"))

## Estimated benefits are the product of estimated probability of interacting with a fire times its estimated benefits: i.e. lambda * mu

facts_act_grouped_06_23 <- mutate(facts_act_grouped_06_23, Benefits = mu*lambda, 
                                  Benefit_Cost = Benefits/TOT_COST)

## Identify observations with high cost/acre

mean_cost <- mean(facts_west$COST_PER_UOM, na.rm = TRUE)
sd_cost <- sd(facts_west$COST_PER_UOM, na.rm = TRUE)
threshold <- mean_cost + 10*sd_cost

## Remove observations with misreported costs from the analysis - i.e. those with cost/acre above 10 SD of the mean cost/acre

facts_act_grouped_06_23_filt <- filter(facts_act_grouped_06_23, COST_PER_UOM <= threshold)

ex <- filter(facts_act_grouped_06_23, ACTIVITY_UNIT_CN == "6150108010602") # one treatment that cost 1 billion USD

## Filter treatments that are above 75 acres in analysis (this makes the BCR smaller but our analysis is based on identifying treatment effects of sufficiently large treatments).

facts_act_grouped_06_23_filt_1 <- filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q_05)

BCR_Matching <- sum(facts_act_grouped_06_23_filt_1$Benefits)/sum(facts_act_grouped_06_23_filt_1$TOT_COST)
N_Treats_Matching <- sum(Treat_Size_Benefits$N_Treats)

Robustness <- read_csv(here("data", "temp", "BCR_Robustness_up.csv"))

Robustness <- rbind(Robustness,c(
  round(BCR_Matching, digits = 2),
  N_Treats_Matching,
  "2017-2023"
))

Robustness <- as.data.frame(Robustness)

rownames(Robustness) <- c("Baseline",  "Alternative Cost Estimate", "No Smoke Exposure Imputation", 
                          "2017-2020 Fires No Imputation", "2017-2020 Fires Imputation", "Yet-to-be Treated Controls Only", "Matched Controls Only")

notes <- "\\parbox[t]{\\textwidth}{\\footnotesize{The following table shows the sensitivity of the ex-ante benefit cost ratio estimates.}}"

# Convert to a single line string
notes_single_line <- gsub("\n", " ", notes)

#### Save Table S10

stargazer(Robustness, 
          title = ("\\label{tab:RobustnessBCR} Sentivitiy of ex-ante benefit cost ratio estimates"), 
          summary = F,
          digits = 0,
          type = "latex", 
          rownames = F,
          colnames = T, 
          out = here("output", "tables", "TableS10.tex"),
          notes = notes_single_line,
          column.sep.width = "-4pt",
          notes.align = "l"
)

toc()
