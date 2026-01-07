##########################################    03_cumulative_effects.R  ##########################################

################# Purpose: From main specification create survival plots 

rm(list=ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr,sf, tmap, magrittr, rnaturalearth, rnaturalearthdata, ggplot2, maps, lwgeom, rgeos, raster, stars, haven, stargazer, quantmod, lubridate, tidyr, ggpubr, 
               rgdal, exactextractr, tictoc, terra, gtools, here, fixest, modelsummary, readr, rdrobust, prism, parallel,tmaptools, 
               OpenStreetMap, maptiles, gifski, geosphere, did, cowplot, didimputation, boot, survival, survminer, patchwork, gridExtra)

# Set Path
here::i_am("code/06_analysis/cumulative_effects.R")


# Load Functions 

source(here("code", "functions","tidy_facts.R"))
source(here("code", "functions","calculate_size_survival.R"))
st_erase = function(x, y) st_difference(x, st_make_valid(st_union(st_combine(y)))) # Taken from here https://r-spatial.github.io/sf/reference/geos_binary_ops.html

## Load in DiD Grids 

Grids_df <- read_csv(here("data", "intermediate", "SpatialDiD_Grids_L24_K05.csv"))
Grids_df <- mutate(Grids_df, Treat_Post = ifelse(Treated_Dir == 1 & time_to_treat >= 0, 1, 0))
Grids_df <- mutate(Grids_df,  dist_treated = ifelse(Treated_Dir == 0, 1000, L_TAU)) # For Sun & Abraham Estimator
Grids_df <- Grids_df %>%
  dplyr::select(
    ID,
    direction_fire_FE,
    distance_bin,
    BURN,
    BURN_LAG,
    time_to_treat,
    everything()
  )  


# Grids_df$LOG_FIRE_INTENSITY[Grids_df$LOG_FIRE_INTENSITY == -Inf] <- NA


#### To do the survival analysis we use an OLS model on the untreated (control & yet-to-be treated) 
##      observations to predict/impute the conditional probability of fire spread for the treated units in 
##      the absence of treatment. We then compare this to the observed survival rates to impute the difference
##      between treated and untreated outcomes.


Grids_df_filt <- filter(Grids_df, BURN_LAG == 1)
Grids_df_filt_1 <- filter(Grids_df_filt, time_to_treat %in% seq(-5,4,1) & Burn_Right_Away == 0)


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


Grids_df_untreated <- filter(Grids_df_filt_1, Treated == 0 & Treated_LAG == 0)
Grids_df_untreated_BURN <- filter(Grids_df_untreated, BURN == 1)

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

# # Dynamically construct the formula
# formula_burn_pred_NAs <- as.formula(paste(
#   "BURN ~ ",
#   paste(controls, collapse = " + "),           # Add all control variables
#   "| distance_bin + FIRE_ID + FuelType_2001 + Aspect_Class"         # Add fixed effects
# ))
# 
# # Dynamically construct the formula
# formula_burn_sev_NAs <- as.formula(paste(
#   "BURN_SEV0 ~ ",
#   paste(controls, collapse = " + "),           # Add all control variables
#   "| distance_bin + FIRE_ID + FuelType_2001 + Aspect_Class"         # Add fixed effects
# ))
# 
# OLS_Pred_Model <- feols(formula_burn_pred,
#                         cluster = ~FIRE_ID,
#                         weights = ~Grid_Acres,
#                         data = Grids_df_untreated)
# 
# OLS_Pred_Model_Burn <- feols(formula_burn_sev,
#                              cluster = ~FIRE_ID,
#                              weights = ~Grid_Acres,
#                              data = Grids_df_untreated_BURN)
# 
# OLS_Pred_Model_Med_High_Burn <- feols(formula_med_high_burn_sev,
#                                       cluster = ~FIRE_ID,
#                                       weights = ~Grid_Acres,
#                                       data = Grids_df_untreated_BURN)
# 
# OLS_Pred_Model_High_Burn <- feols(formula_high_burn_sev,
#                                   cluster = ~FIRE_ID,
#                                   weights = ~Grid_Acres,
#                                   data = Grids_df_untreated_BURN)
# 
# 
# OLS_Pred_Model_NAs <- feols(formula_burn_pred_NAs,
#                         cluster = ~FIRE_ID,
#                         weights = ~Grid_Acres,
#                         data = Grids_df_untreated)
# 
# OLS_Pred_Model_Burn_NAs <- feols(formula_burn_sev_NAs,
#                              cluster = ~FIRE_ID,
#                              weights = ~Grid_Acres,
#                              data = Grids_df_untreated_BURN)



# #### Aside what is R^2 with and without MTT outputs - addressing reviewer 1's concern ####
# 
# Grids_df_untreated$LOG_FIRE_INTENSITY[Grids_df_untreated$LOG_FIRE_INTENSITY == -Inf] <- NA
# Grids_df_untreated_BURN <- filter(Grids_df_untreated, BURN == 1)
# 
# 
# mod1 <- feols(scale(BURN) ~ scale(DELTA_MTT) | direction_fire_FE + distance_bin,
#       cluster = ~FIRE_ID,
#       weights = ~Grid_Acres,
#       data = Grids_df_untreated)
# 
# summary(mod1)
# 
# mod2 <- feols(scale(BURN) ~ scale(LOG_FIRE_INTENSITY) | direction_fire_FE + distance_bin,
#                                cluster = ~FIRE_ID,
#                                weights = ~Grid_Acres,
#                                data = Grids_df_untreated)
# 
# summary(mod2)
# 
# mod3 <- feols(scale(BURN_SEV_MED_HIGH_PCT) ~ scale(DELTA_MTT) | direction_fire_FE + distance_bin,
#       cluster = ~FIRE_ID,
#       weights = ~Grid_Acres,
#       data = Grids_df_untreated)
# 
# summary(mod3)
# 
# mod4 <- feols(scale(BURN_SEV_MED_HIGH_PCT) ~ scale(LOG_FIRE_INTENSITY) | direction_fire_FE + distance_bin,
#       cluster = ~FIRE_ID,
#       weights = ~Grid_Acres,
#       data = Grids_df_untreated)
# 
# summary(mod4)
# 
# ####           Return                     ####




time_window <- seq(-1,9,1)

treated_data <- Grids_df %>% filter(Treated_Dir == 1 & Burn_Right_Away == 0) %>% filter(time_to_treat %in% time_window)

treated_data$pred_without_treatment <- predict(OLS_Pred_Model, newdata = treated_data, type = "response")
treated_data$pred_without_treatment_BS <- predict(OLS_Pred_Model_Burn, newdata = treated_data, type = "response")
treated_data$pred_without_treatment_MED_HIGH_BS <- predict(OLS_Pred_Model_Med_High_Burn, newdata = treated_data, type = "response")
treated_data$pred_without_treatment_HIGH_BS <- predict(OLS_Pred_Model_High_Burn, newdata = treated_data, type = "response")


sum(is.na(treated_data$pred_without_treatment))
sum(is.na(treated_data$pred_without_treatment_BS))
sum(is.na(treated_data$pred_without_treatment_MED_HIGH_BS))
sum(is.na(treated_data$pred_without_treatment_HIGH_BS))


# nas <- dplyr::select(treated_data, direction_fire_FE, BURN, time_to_treat, dist_treated, distance_bin, pred_without_treatment, pred_without_treatment_BS) %>% arrange(direction_fire_FE) %>%
#   filter(is.na(pred_without_treatment_BS) == T)
# 
# Grids_df_untreated[Grids_df_untreated$direction_fire_FE == "3-CA4172612124320190728",]
# ex <- Grids_df[Grids_df$direction_fire_FE == "3-CA4172612124320190728",]


# 1. Get alternate predictions from the Fire-FE model for all rows
# treated_data <- treated_data %>%
#   mutate(pred_without_treatment_alt = predict(
#     OLS_Pred_Model_NAs,
#     newdata = cur_data_all(),   # dplyr >= 1.1.0
#     type   = "response"
#   ))
# 
# # 2. Replace NA's in the main prediction with the alt prediction
# treated_data <- treated_data %>%
#   mutate(pred_without_treatment = if_else(
#     is.na(pred_without_treatment),
#     pred_without_treatment_alt,
#     pred_without_treatment
#   )) %>%
#   dplyr::select(-pred_without_treatment_alt)  # optional cleanup

treated_data <- filter(treated_data, is.na(pred_without_treatment) == F)

sum(is.na(treated_data$pred_without_treatment))

ex <- dplyr::select(treated_data, direction_fire_FE, BURN, time_to_treat, dist_treated, distance_bin, pred_without_treatment) %>% arrange(direction_fire_FE)

treated_data[is.na(treated_data$pred_without_treatment), "pred_without_treatment"] <- 0

## Option 1 - for predictions > 1 set to 1, < 0, set to 0

# treated_data <- mutate(treated_data, pred_without_treatment = case_when(pred_without_treatment > 1 ~ 1,
#                                                                      pred_without_treatment < 0 ~ 0,
#                                                                      pred_without_treatment >= 0 & pred_without_treatment <= 1 ~ pred_without_treatment))


treated_data[is.na(treated_data$pred_without_treatment_BS), "pred_without_treatment_BS"] <- 0
treated_data[is.na(treated_data$pred_without_treatment_MED_HIGH_BS), "pred_without_treatment_MED_HIGH_BS"] <- 0
treated_data[is.na(treated_data$pred_without_treatment_HIGH_BS), "pred_without_treatment_HIGH_BS"] <- 0

# treated_data <- mutate(treated_data, pred_without_treatment_HIGH_BS = case_when(pred_without_treatment_HIGH_BS > 1 ~ 1,
#                                                                                 pred_without_treatment_HIGH_BS < 0 ~ 0,
#                                                                                 pred_without_treatment_HIGH_BS >= 0 & pred_without_treatment_HIGH_BS <= 1 ~ pred_without_treatment_HIGH_BS))

treated_data <- mutate(treated_data, 
                       pred_without_treatment_MED_HIGH_BS = case_when(pred_without_treatment_MED_HIGH_BS > 1 ~ 1,
                                                                      pred_without_treatment_MED_HIGH_BS < 0 ~ 0,
                                                                      pred_without_treatment_MED_HIGH_BS >= 0 & pred_without_treatment_MED_HIGH_BS <= 1 ~ pred_without_treatment_MED_HIGH_BS),
                       pred_without_treatment_HIGH_BS = case_when(pred_without_treatment_HIGH_BS > 1 ~ 1,
                                                                  pred_without_treatment_HIGH_BS < 0 ~ 0,
                                                                  pred_without_treatment_HIGH_BS >= 0 & pred_without_treatment_HIGH_BS <= 1 ~ pred_without_treatment_HIGH_BS))

treated_data <- mutate(treated_data, 
                       pred_with_treatment = BURN,
                       pred_with_treatment_BS = BURN_SEV0,
                       pred_with_treatment_MED_HIGH_BS = BURN_SEV_MED_HIGH_PCT,
                       pred_with_treatment_HIGH_BS = BURN_SEV_HIGH_PCT) %>%
  arrange(direction_fire_FE, time_to_treat)


# Create 1-6 lags for each variable
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


N_obs_treat_time <- treated_data %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n())

N_Directions <- as.numeric(N_obs_treat_time[1,2])

survivor_plot <- treated_data %>%
  group_by(time_to_treat) %>%
  summarise(Mean_Survival_Treat = (sum(survival_with_treatment)/N_Directions),
            Mean_Survival_No_Treat = (sum(survival_without_treatment)/N_Directions),
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

survivor_plot <- mutate(survivor_plot, Pct_BS_Change = (Mean_BS_No_Treat - Mean_BS_Treat)/Mean_BS_No_Treat,
                        Pct_MED_HIGH_BS_Change = (Pct_Burned_MS_HS_Treat - Pct_Burned_MS_HS_No_Treat)/Pct_Burned_MS_HS_No_Treat,
                        Pct_HIGH_BS_Change = (Pct_Burned_HS_Treat - Pct_Burned_HS_No_Treat)/Pct_Burned_HS_No_Treat
)

survivor_plot$distance <- seq(0, 2.5, 0.5)

interval_labels <- c(
  "0", "0.5", "1",
  "1.5", "2", "2.5"
)

# Plotting the time series for predicted survival rates
survival_plot <- ggplot(survivor_plot, aes(x = distance)) +
  # Dashed vertical line at -0.5
  geom_line(aes(y = Mean_Survival_Treat, color = "With Treatment (Observed)"), size = 1.2) +
  geom_line(aes(y = Mean_Survival_No_Treat, color = "Without Treatment (Counterfactual)"), size = 1.2) +
  labs(
    title = "Burn Probability by Distance From Treatment",
    x = "Distance from treatment interaction (km)",
    y = "Mean Predicted Burn Probability",
    color = "Treatment Status"
  ) +
  scale_color_manual(values = c("With Treatment (Observed)" = "#0072B2", "Without Treatment (Counterfactual)" = "#D55E00")) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text  = element_text(size = 7),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size  = 6,      # smaller than the y-axis labels
      lineheight = 0.9
    ),
    plot.title   = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                margin = margin(b = 2)),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(3, 3, 3, 3)
  )

survival_plot

# Plotting the time series for predicted survival rates
survival_plot_BS <- ggplot(survivor_plot, aes(x = time_to_treat)) +
  # Dashed vertical line at -0.5
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_line(aes(y = Mean_BS_Treat, color = "With Treatment"), size = 1.2) +
  geom_line(aes(y = Mean_BS_No_Treat, color = "Without Treatment"), size = 1.2) +
  labs(
    title = "Predicted Burn Severity by Distance From Treatment",
    x = "Distance To Treatment (0.5km)",
    y = "Mean Predicted Burn Severity",
    color = "Treatment Status"
  ) +
  scale_color_manual(values = c("With Treatment" = "#0072B2", "Without Treatment" = "#D55E00")) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),  # Remove minor gridlines for a cleaner look
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add axes
    axis.line = element_blank(),  # Avoid duplicate axis lines
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

# Plotting the time series for predicted survival rates
survival_plot_MED_HIGH_BS <- ggplot(survivor_plot, aes(x = time_to_treat)) +
  # Dashed vertical line at -0.5
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_line(aes(y = Pct_Burned_MS_HS_Treat, color = "With Treatment"), size = 1.2) +
  geom_line(aes(y = Pct_Burned_MS_HS_No_Treat, color = "Without Treatment"), size = 1.2) +
  labs(
    title = "Predicted Share of Med-High Burn Severity by Distance From Treatment",
    x = "Distance To Treatment (0.5km)",
    y = "Mean Predicted Share of Med-High Burn Severity",
    color = "Treatment Status"
  ) +
  scale_color_manual(values = c("With Treatment" = "#0072B2", "Without Treatment" = "#D55E00")) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),  # Remove minor gridlines for a cleaner look
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add axes
    axis.line = element_blank(),  # Avoid duplicate axis lines
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

# survival_plot_HIGH_BS <- ggplot(survivor_plot, aes(x = time_to_treat)) +
#   # Dashed vertical line at -0.5
#   geom_vline(xintercept = -0.5, linetype = "dashed", color = "black", linewidth = 0.8) +
#   geom_line(aes(y = Pct_Burned_HS_Treat, color = "With Treatment"), size = 1.2) +
#   geom_line(aes(y = Pct_Burned_HS_No_Treat, color = "Without Treatment"), size = 1.2) +
#   labs(
#     title = "Predicted Share of High Burn Severity by Distance From Treatment",
#     x = "Distance To Treatment (0.5km)",
#     y = "Mean Predicted Share of High Burn Severity",
#     color = "Treatment Status"
#   ) +
#   scale_color_manual(values = c("With Treatment" = "#0072B2", "Without Treatment" = "#D55E00")) +
#   theme_minimal(base_size = 14) +
#   theme(
#     panel.grid = element_blank(),  # Remove minor gridlines for a cleaner look
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add axes
#     axis.line = element_blank(),  # Avoid duplicate axis lines
#     axis.title = element_text(size = 16, face = "bold"),
#     axis.text = element_text(size = 16),
#     plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

bootstrap_survivor_LPM <- function(data, indices) {
  # Resample clusters (FIRE_ID)
  resampled_fires <- unique(data$FIRE_ID)[indices]
  
  Resample_Grids <- filter(data, FIRE_ID %in% resampled_fires)
  
  Grids_df_filt <- filter(Resample_Grids, BURN_LAG == 1)
  Grids_df_filt_1 <- filter(Grids_df_filt, time_to_treat %in% seq(-5,4,1) & Burn_Right_Away == 0)
  
  Grids_df_untreated <- filter(Grids_df_filt_1, Treated == 0 & Treated_LAG == 0)
  Grids_df_untreated_BURN <- filter(Grids_df_untreated, BURN == 1)
  
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
  
  
  treated_data <- Resample_Grids %>% filter(Treated_Dir == 1 & Burn_Right_Away == 0) %>% filter(time_to_treat %in% seq(-1,4,1))
  
  treated_data$pred_without_treatment <- predict(OLS_Pred_Model, newdata = treated_data, type = "response")
  treated_data$pred_without_treatment_BS <- predict(OLS_Pred_Model_Burn, newdata = treated_data, type = "response")
  treated_data$pred_without_treatment_MED_HIGH_BS <- predict(OLS_Pred_Model_Med_High_Burn, newdata = treated_data, type = "response")
  treated_data$pred_without_treatment_HIGH_BS <- predict(OLS_Pred_Model_High_Burn, newdata = treated_data, type = "response")
  
  treated_data <- filter(treated_data, is.na(pred_without_treatment) == F)
  
  
  treated_data[is.na(treated_data$pred_without_treatment), "pred_without_treatment"] <- 0
  treated_data[is.na(treated_data$pred_without_treatment_BS), "pred_without_treatment_BS"] <- 0
  treated_data[is.na(treated_data$pred_without_treatment_MED_HIGH_BS), "pred_without_treatment_MED_HIGH_BS"] <- 0
  treated_data[is.na(treated_data$pred_without_treatment_HIGH_BS), "pred_without_treatment_HIGH_BS"] <- 0
  
  # treated_data <- mutate(treated_data, pred_without_treatment = case_when(pred_without_treatment > 1 ~ 1,
  #                                                                         pred_without_treatment < 0 ~ 0,
  #                                                                         pred_without_treatment >= 0 & pred_without_treatment <= 1 ~ pred_without_treatment))
  
  treated_data <- mutate(treated_data, 
                         pred_without_treatment_MED_HIGH_BS = case_when(pred_without_treatment_MED_HIGH_BS > 1 ~ 1,
                                                                        pred_without_treatment_MED_HIGH_BS < 0 ~ 0,
                                                                        pred_without_treatment_MED_HIGH_BS >= 0 & pred_without_treatment_MED_HIGH_BS <= 1 ~ pred_without_treatment_MED_HIGH_BS),
                         pred_without_treatment_HIGH_BS = case_when(pred_without_treatment_HIGH_BS > 1 ~ 1,
                                                                    pred_without_treatment_HIGH_BS < 0 ~ 0,
                                                                    pred_without_treatment_HIGH_BS >= 0 & pred_without_treatment_HIGH_BS <= 1 ~ pred_without_treatment_HIGH_BS))
  
  treated_data <- mutate(treated_data, 
                         pred_with_treatment = BURN,
                         pred_with_treatment_BS = BURN_SEV0,
                         pred_with_treatment_MED_HIGH_BS = BURN_SEV_MED_HIGH_PCT,
                         pred_with_treatment_HIGH_BS = BURN_SEV_HIGH_PCT) %>%
    arrange(direction_fire_FE, time_to_treat)
  
  
  # Create 1-6 lags for each variable
  treated_data <- treated_data %>%
    arrange(direction_fire_FE, time_to_treat) %>% # Ensure data is ordered by ID and time
    group_by(direction_fire_FE) %>%              # Group by ID
    mutate(p1 = lag(pred_with_treatment, 1),
           p2 = lag(pred_with_treatment, 2),
           p3 = lag(pred_with_treatment, 3),
           p4 = lag(pred_with_treatment, 4),
           p10 = lag(pred_without_treatment, 1),
           p20 = lag(pred_without_treatment, 2),
           p30 = lag(pred_without_treatment, 3),
           p40 = lag(pred_without_treatment, 4)) %>% # Create lagged columns
    ungroup()                                      # Remove grouping
  
  treated_data <- mutate(treated_data, 
                         survival_with_treatment = case_when(time_to_treat == -1 ~ 1,
                                                             time_to_treat == 0 ~ pred_with_treatment,
                                                             time_to_treat == 1 ~ pred_with_treatment*p1, 
                                                             time_to_treat == 2 ~ pred_with_treatment*p1*p2,
                                                             time_to_treat == 3 ~ pred_with_treatment*p1*p2*p3,
                                                             time_to_treat == 4 ~ pred_with_treatment*p1*p2*p3*p4),
                         survival_without_treatment = case_when(time_to_treat == -1 ~ 1,
                                                                time_to_treat == 0 ~ pred_without_treatment,
                                                                time_to_treat == 1 ~ pred_without_treatment*p10, 
                                                                time_to_treat == 2 ~ pred_without_treatment*p10*p20,
                                                                time_to_treat == 3 ~ pred_without_treatment*p10*p20*p30,
                                                                time_to_treat == 4 ~ pred_without_treatment*p10*p20*p30*p40)
  )
  
  
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
  
  N_obs_treat_time <- treated_data %>%
    group_by(time_to_treat) %>%
    summarise(N_obs = n())
  
  
  N_Directions <- as.numeric(N_obs_treat_time[1,2])
  
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
              Diff_MED_HIGH_BS = Pct_Burned_MS_HS_Treat - Pct_Burned_MS_HS_No_Treat,
              Diff_HIGH_BS = Pct_Burned_HS_Treat - Pct_Burned_HS_No_Treat,
              Acres_Burned_MS_HS_Treat = sum(MED_HIGH_BS_with_treatment*Grid_Acres),
              Acres_Burned_MS_HS_No_Treat = sum(MED_HIGH_BS_without_treatment*Grid_Acres),
              Acres_Burned_HS_Treat = sum(HIGH_BS_with_treatment*Grid_Acres),
              Acres_Burned_HS_No_Treat = sum(HIGH_BS_without_treatment*Grid_Acres)) %>%
    filter(time_to_treat %in% seq(0,5,1))
  
  leads_Treat <- as.numeric(survivor_plot$Mean_Survival_Treat)
  leads_No_Treat <- as.numeric(survivor_plot$Mean_Survival_No_Treat)
  
  if (type == "Ordinal"){
    leads_Treat_BS <- as.numeric(survivor_plot$Mean_BS_Treat)
    leads_No_Treat_BS <- as.numeric(survivor_plot$Mean_BS_No_Treat)
    diff_BS <- as.numeric(survivor_plot$Diff_BS)
  }
  
  if (type == "Med-High"){
    leads_Treat_BS <- as.numeric(survivor_plot$Pct_Burned_MS_HS_Treat)
    leads_No_Treat_BS <- as.numeric(survivor_plot$Pct_Burned_MS_HS_No_Treat)
    diff_BS <- as.numeric(survivor_plot$Diff_MED_HIGH_BS)
  }
  
  if (type == "High"){
    leads_Treat_BS <- as.numeric(survivor_plot$Pct_Burned_HS_Treat)
    leads_No_Treat_BS <- as.numeric(survivor_plot$Pct_Burned_HS_No_Treat)
    diff_BS <- as.numeric(survivor_plot$Diff_HIGH_BS)
  }
  
  diff_survival <- as.numeric(survivor_plot$Diff_Survival)
  
  return(c(leads_Treat, leads_No_Treat, leads_Treat_BS, leads_No_Treat_BS, diff_survival, diff_BS))
}

set.seed(123) # For reproducibility
type <- "Med-High"
bootstrap_results <- boot(data = Grids_df, 
                          statistic = bootstrap_survivor_LPM, 
                          parallel = "multicore", 
                          ncpus = 8,
                          R = 1000)


# Number of leads for treated and untreated
num_leads <- 5

# Extract bootstrap replicates
bootstrap_replicates <- bootstrap_results$t

# Initialize matrices to store confidence intervals
treated_ci <- matrix(NA, nrow = num_leads, ncol = 2)
untreated_ci <- matrix(NA, nrow = num_leads, ncol = 2)
diff_ci <- matrix(NA, nrow = num_leads, ncol = 2)



treated_BS_ci <- matrix(NA, nrow = num_leads, ncol = 2)
untreated_BS_ci <- matrix(NA, nrow = num_leads, ncol = 2)
diff_BS_ci <- matrix(NA, nrow = num_leads, ncol = 2)

# Confidence level (e.g., 95%)
alpha <- 0.05

# Compute confidence intervals for treated and untreated leads
for (i in 1:num_leads) {
  # Treated
  treated_ci[i, ] <- quantile(bootstrap_replicates[, i], probs = c(alpha / 2, 1 - alpha / 2))
  
  # Untreated
  untreated_ci[i, ] <- quantile(bootstrap_replicates[, i + num_leads], probs = c(alpha / 2, 1 - alpha / 2))
  
  # Treated
  treated_BS_ci[i, ] <- quantile(bootstrap_replicates[, 10 + i], probs = c(alpha / 2, 1 - alpha / 2))
  
  # Untreated
  untreated_BS_ci[i, ] <- quantile(bootstrap_replicates[, 10 + i + num_leads], probs = c(alpha / 2, 1 - alpha / 2))
  
  diff_ci[i, ] <- quantile(bootstrap_replicates[, 20 + i], probs = c(alpha / 2, 1 - alpha / 2))
  diff_BS_ci[i, ] <- quantile(bootstrap_replicates[, 20 + i + num_leads], probs = c(alpha / 2, 1 - alpha / 2))
  
  
}

# Combine results into a data frame for easy visualization
ci_results <- data.frame(
  Lead = 0:(num_leads - 1),
  Treated_Mean = bootstrap_results$t0[1:num_leads],
  Treated_Lower = treated_ci[, 1],
  Treated_Upper = treated_ci[, 2],
  Untreated_Mean = bootstrap_results$t0[(num_leads + 1):(2 * num_leads)],
  Untreated_Lower = untreated_ci[, 1],
  Untreated_Upper = untreated_ci[, 2],
  
  Treated_BS_Mean = bootstrap_results$t0[11:15],
  Treated_BS_Lower = treated_BS_ci[, 1],
  Treated_BS_Upper = treated_BS_ci[, 2],
  Untreated_BS_Mean = bootstrap_results$t0[16:20],
  Untreated_BS_Lower = untreated_BS_ci[, 1],
  Untreated_BS_Upper = untreated_BS_ci[, 2],
  
  Diff_Mean = bootstrap_results$t0[21:25],
  Diff_Lower = diff_ci[, 1],
  Diff_Upper = diff_ci[, 2],
  
  Diff_BS_Mean = bootstrap_results$t0[26:30],
  Diff_BS_Lower = diff_BS_ci[, 1],
  Diff_BS_Upper = diff_BS_ci[, 2]
)

# Print results
print(ci_results)


survivor_plot_boot <- merge(survivor_plot, ci_results, by.x = "time_to_treat", by.y = "Lead", all.x = T)

if (type == "Ordinal"){
  starting_BS <- filter(survivor_plot_boot, time_to_treat == -1)$Mean_BS_Treat
}

if (type == "Med-High"){
  starting_BS <- filter(survivor_plot_boot, time_to_treat == -1)$Pct_Burned_MS_HS_Treat
}

if (type == "High"){
  starting_BS <- filter(survivor_plot_boot, time_to_treat == -1)$Pct_Burned_HS_Treat
}

survivor_plot_boot[1, c("Treated_Mean", "Treated_Lower", "Treated_Upper", 
                        "Untreated_Mean", "Untreated_Lower", "Untreated_Upper")] <- rep(1, 6)

survivor_plot_boot[1, c("Diff_Mean", "Diff_Lower", "Diff_Upper", 
                        "Diff_BS_Mean", "Diff_BS_Lower", "Diff_BS_Upper")] <- rep(0, 6)


survivor_plot_boot[1, c("Treated_BS_Mean", "Treated_BS_Lower", "Treated_BS_Upper", 
                        "Untreated_BS_Mean", "Untreated_BS_Lower", "Untreated_BS_Upper")] <- rep(starting_BS, 6)

# survivor_plot_boot <- filter(survivor_plot_boot, time_to_treat >= 0)
# survivor_plot_boot$distance <- seq(0, 2.5, 0.5)



survival_plot_boot <- ggplot(survivor_plot_boot, aes(x = distance)) +
  # Confidence intervals as ribbons
  geom_ribbon(aes(ymin = Treated_Lower*100, ymax = Treated_Upper*100, fill = "With Treatment (Observed)"), alpha = 0.2) +
  geom_ribbon(aes(ymin = Untreated_Lower*100, ymax = Untreated_Upper*100, fill = "Without Treatment (Counterfactual)"), alpha = 0.2) +
  # Mean survival lines
  geom_line(aes(y = Mean_Survival_Treat*100, color = "With Treatment (Observed)"), size = 1.2) +
  geom_line(aes(y = Mean_Survival_No_Treat*100, color = "Without Treatment (Counterfactual)"), size = 1.2) +
  
  # Labels and titles
  labs(
    title = "Cumulative Burn Probability",
    x = "Distance from treatment interaction (km)",
    y = "Predicted burn probability (%)",
    color = "Treatment Status",
    fill = "Treatment Status"
  ) +
  
  # Custom colors for lines and ribbons
  scale_color_manual(values = c("With Treatment (Observed)" = "#0072B2", "Without Treatment (Counterfactual)" = "#D55E00")) +
  scale_fill_manual(values = c("With Treatment (Observed)" = "#0072B2", "Without Treatment (Counterfactual)" = "#D55E00")) +
  
  # Theme adjustments
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text  = element_text(size = 7),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size  = 6,      # smaller than the y-axis labels
      lineheight = 0.9
    ),
    plot.title   = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                margin = margin(b = 2)),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(3, 3, 3, 3)
  )


survival_diff_plot_boot <- ggplot(survivor_plot_boot, aes(x = distance)) +
  # Confidence intervals as ribbons
  geom_ribbon(aes(ymin = Diff_Lower*100, ymax = Diff_Upper*100, fill = "Observed - Counterfactual (%)"), alpha = 0.2) +
  # Mean survival lines
  geom_line(aes(y = Diff_Mean*100, color = "Observed - Counterfactual (%)"), size = 1.2) +
  # Dotted horizontal line at zero
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  # Labels and titles
  labs(
    title = "Cumulative Effects on Burn Probability",
    x = "Distance from treatment interaction (km)",
    y = "Reduction in burn probability (% points)",
    color = "Difference",
    fill = "Difference"
  ) +
  
  # Custom colors for lines and ribbons
  scale_color_manual(values = c("Observed - Counterfactual (%)" = "#0072B2")) +
  
  # Theme adjustments
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text  = element_text(size = 7),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size  = 6,      # smaller than the y-axis labels
      lineheight = 0.9
    ),
    plot.title   = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                margin = margin(b = 2)),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(3, 3, 3, 3)
  )

if (type == "Med-High"){
  # Create survivor_BS_plot with no legend
  survival_MH_BS_plot_boot <- ggplot(survivor_plot_boot, aes(x = distance)) +
    # Confidence intervals as ribbons
    geom_ribbon(aes(ymin = 100*Treated_BS_Lower, ymax = Treated_BS_Upper*100, fill = "With Treatment (Observed)"), alpha = 0.2) +
    geom_ribbon(aes(ymin = 100*Untreated_BS_Lower, ymax = Untreated_BS_Upper*100, fill = "Without Treatment (Counterfactual)"), alpha = 0.2) +
    # Mean survival lines
    geom_line(aes(y = Pct_Burned_MS_HS_Treat*100, color = "With Treatment (Observed)"), size = 1.2) +
    geom_line(aes(y = Pct_Burned_MS_HS_No_Treat*100, color = "Without Treatment (Counterfactual)"), size = 1.2) +
    
    # Labels and titles
    labs(
      title = "Cumulative Moderate-High Severity Burn",
      x = "Distance from treatment interaction (km)",
      y = "Predicted moderate–high severity burn share (%)",
      color = "Treatment Status",
      fill = "Treatment Status"
    ) +
    
    # Custom colors for lines and ribbons
    scale_color_manual(values = c("With Treatment (Observed)" = "#0072B2", "Without Treatment (Counterfactual)" = "#D55E00")) +
    scale_fill_manual(values = c("With Treatment (Observed)" = "#0072B2", "Without Treatment (Counterfactual)" = "#D55E00")) +
    
    # Theme adjustments
    theme_minimal(base_family = "Helvetica") +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.title = element_text(size = 8, face = "bold"),
      axis.text  = element_text(size = 7),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size  = 6,      # smaller than the y-axis labels
        lineheight = 0.9
      ),
      plot.title   = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                  margin = margin(b = 2)),
      axis.ticks = element_line(linewidth = 0.3),
      plot.margin = margin(3, 3, 3, 3)
    )
  
}

if (type == "Med-High"){
  survival_MH_BS_diff_plot_boot <- ggplot(survivor_plot_boot, aes(x = distance)) +
    # Confidence intervals as ribbons
    geom_ribbon(aes(ymin = Diff_BS_Lower*100, ymax = Diff_BS_Upper*100, fill = "Observed - Counterfactual (%)"), alpha = 0.2) +
    # Mean survival lines
    geom_line(aes(y = Diff_BS_Mean*100, color = "Observed - Counterfactual (%)"), size = 1.2) +
    # Dotted horizontal line at zero
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    # Labels and titles
    labs(
      title = "Cumulative Effects on Moderate-High Severity Burn",
      x = "Distance from treatment interaction (km)",
      y = "Reduction in moderate–high severity burn share (% points)",
      color = "Difference",
      fill = "Difference"
    ) +
    
    # Custom colors for lines and ribbons
    scale_color_manual(values = c("Observed - Counterfactual (%)" = "#0072B2")) +
    
    # Theme adjustments
    theme_minimal(base_family = "Helvetica") +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.title = element_text(size = 8, face = "bold"),
      axis.text  = element_text(size = 7),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size  = 6,      # smaller than the y-axis labels
        lineheight = 0.9
      ),
      plot.title   = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                  margin = margin(b = 2)),
      axis.ticks = element_line(linewidth = 0.3),
      plot.margin = margin(3, 3, 3, 3)
    )}



# Combine the two plots and ensure only one legend
final_layout <- (survival_plot_boot | survival_MH_BS_plot_boot) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(
      family = "Helvetica",
      face   = "bold",
      size   = 10
    ),
    legend.title = element_text(
      family = "Helvetica",
      face   = "bold",
      size   = 7
    ),
    legend.text = element_text(
      family = "Helvetica",
      size   = 7
    ),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "bottom"
  )


# plot1 <- readRDS(here("output", "event_plotBURN.rds"))
# plot2 <- readRDS(here("output", "event_plotBURNSEV.rds"))

ggsave(here("output", "figures", "Figure4.pdf"), final_layout,
       width  = 7.24,          # Science 3-column width
       height = 4.5,    # ~4.5 in
       units  = "in",
       dpi    = 300)

# Combine the two plots and ensure only one legend
final_layout <- (survival_diff_plot_boot | survival_MH_BS_diff_plot_boot) +
  plot_layout(guides = "collect") +  # <- collect shared legend
  plot_annotation(tag_levels = "A")  &
  theme(
    plot.tag = element_text(
      family = "Helvetica",
      face   = "bold",
      size   = 10
    ),
    legend.title = element_text(
      family = "Helvetica",
      face   = "bold",
      size   = 7
    ),
    legend.text = element_text(
      family = "Helvetica",
      size   = 7
    ),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "bottom"
  )

ggsave(here("output", "figures", "Survival_Difference.pdf"), final_layout, 
       width  = 7.24,          # Science 3-column width
       height = 4.5,    # ~4.5 in
       units  = "in",
       dpi    = 300)


###### Ordinal Burn Severity

set.seed(123) # For reproducibility
type <- "Ordinal"
bootstrap_results <- boot(data = Grids_df, 
                          statistic = bootstrap_survivor_LPM, 
                          parallel = "multicore", 
                          ncpus = 8,
                          R = 1000)


# Number of leads for treated and untreated
num_leads <- 5

# Extract bootstrap replicates
bootstrap_replicates <- bootstrap_results$t

# Initialize matrices to store confidence intervals
treated_ci <- matrix(NA, nrow = num_leads, ncol = 2)
untreated_ci <- matrix(NA, nrow = num_leads, ncol = 2)
diff_ci <- matrix(NA, nrow = num_leads, ncol = 2)



treated_BS_ci <- matrix(NA, nrow = num_leads, ncol = 2)
untreated_BS_ci <- matrix(NA, nrow = num_leads, ncol = 2)
diff_BS_ci <- matrix(NA, nrow = num_leads, ncol = 2)

# Confidence level (e.g., 95%)
alpha <- 0.05

# Compute confidence intervals for treated and untreated leads
for (i in 1:num_leads) {
  # Treated
  treated_ci[i, ] <- quantile(bootstrap_replicates[, i], probs = c(alpha / 2, 1 - alpha / 2))
  
  # Untreated
  untreated_ci[i, ] <- quantile(bootstrap_replicates[, i + num_leads], probs = c(alpha / 2, 1 - alpha / 2))
  
  # Treated
  treated_BS_ci[i, ] <- quantile(bootstrap_replicates[, 10 + i], probs = c(alpha / 2, 1 - alpha / 2))
  
  # Untreated
  untreated_BS_ci[i, ] <- quantile(bootstrap_replicates[, 10 + i + num_leads], probs = c(alpha / 2, 1 - alpha / 2))
  
  diff_ci[i, ] <- quantile(bootstrap_replicates[, 20 + i], probs = c(alpha / 2, 1 - alpha / 2))
  diff_BS_ci[i, ] <- quantile(bootstrap_replicates[, 20 + i + num_leads], probs = c(alpha / 2, 1 - alpha / 2))
  
  
}

# Combine results into a data frame for easy visualization
ci_results <- data.frame(
  Lead = 0:(num_leads - 1),
  Treated_Mean = bootstrap_results$t0[1:num_leads],
  Treated_Lower = treated_ci[, 1],
  Treated_Upper = treated_ci[, 2],
  Untreated_Mean = bootstrap_results$t0[(num_leads + 1):(2 * num_leads)],
  Untreated_Lower = untreated_ci[, 1],
  Untreated_Upper = untreated_ci[, 2],
  
  Treated_BS_Mean = bootstrap_results$t0[11:15],
  Treated_BS_Lower = treated_BS_ci[, 1],
  Treated_BS_Upper = treated_BS_ci[, 2],
  Untreated_BS_Mean = bootstrap_results$t0[16:20],
  Untreated_BS_Lower = untreated_BS_ci[, 1],
  Untreated_BS_Upper = untreated_BS_ci[, 2],
  
  Diff_Mean = bootstrap_results$t0[21:25],
  Diff_Lower = diff_ci[, 1],
  Diff_Upper = diff_ci[, 2],
  
  Diff_BS_Mean = bootstrap_results$t0[26:30],
  Diff_BS_Lower = diff_BS_ci[, 1],
  Diff_BS_Upper = diff_BS_ci[, 2]
)

# Print results
print(ci_results)


survivor_plot_boot <- merge(survivor_plot, ci_results, by.x = "time_to_treat", by.y = "Lead", all.x = T)

if (type == "Ordinal"){
  starting_BS <- filter(survivor_plot_boot, time_to_treat == -1)$Mean_BS_Treat
}

survivor_plot_boot[1, c("Treated_Mean", "Treated_Lower", "Treated_Upper", 
                        "Untreated_Mean", "Untreated_Lower", "Untreated_Upper")] <- rep(1, 6)

survivor_plot_boot[1, c("Diff_Mean", "Diff_Lower", "Diff_Upper", 
                        "Diff_BS_Mean", "Diff_BS_Lower", "Diff_BS_Upper")] <- rep(0, 6)


survivor_plot_boot[1, c("Treated_BS_Mean", "Treated_BS_Lower", "Treated_BS_Upper", 
                        "Untreated_BS_Mean", "Untreated_BS_Lower", "Untreated_BS_Upper")] <- rep(starting_BS, 6)

if (type == "Ordinal"){
  # Create survivor_BS_plot with no legend
  survival_BS_plot_boot <- ggplot(survivor_plot_boot, aes(x = distance)) +
    # Confidence intervals as ribbons
    geom_ribbon(aes(ymin = Treated_BS_Lower, ymax = Treated_BS_Upper, fill = "With Treatment (Observed)"), alpha = 0.2) +
    geom_ribbon(aes(ymin = Untreated_BS_Lower, ymax = Untreated_BS_Upper, fill = "Without Treatment (Counterfactual)"), alpha = 0.2) +
    # Mean survival lines
    geom_line(aes(y = Mean_BS_Treat, color = "With Treatment (Observed)"), size = 1.2) +
    geom_line(aes(y = Mean_BS_No_Treat, color = "Without Treatment (Counterfactual)"), size = 1.2) +
    
    # Labels and titles
    labs(
      title = "Cumulative Burn Severity",
      x = "Distance from treatment interaction (km)",
      y = "Predicted burn severity",
      color = "Treatment Status",
      fill = "Treatment Status"
    ) +
    
    # Custom colors for lines and ribbons
    scale_color_manual(values = c("With Treatment (Observed)" = "#0072B2", "Without Treatment (Counterfactual)" = "#D55E00")) +
    scale_fill_manual(values = c("With Treatment (Observed)" = "#0072B2", "Without Treatment (Counterfactual)" = "#D55E00")) +
    
    # Theme adjustments
    theme_minimal(base_family = "Helvetica") +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.title = element_text(size = 8, face = "bold"),
      axis.text  = element_text(size = 7),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size  = 6,      # smaller than the y-axis labels
        lineheight = 0.9
      ),
      plot.title   = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                  margin = margin(b = 2)),
      axis.ticks = element_line(linewidth = 0.3),
      plot.margin = margin(3, 3, 3, 3)
    )
}


###### High Burn Severity

set.seed(123) # For reproducibility
type <- "High"
bootstrap_results <- boot(data = Grids_df, 
                          statistic = bootstrap_survivor_LPM, 
                          parallel = "multicore", 
                          ncpus = 8,
                          R = 1000)


# Number of leads for treated and untreated
num_leads <- 5

# Extract bootstrap replicates
bootstrap_replicates <- bootstrap_results$t

# Initialize matrices to store confidence intervals
treated_ci <- matrix(NA, nrow = num_leads, ncol = 2)
untreated_ci <- matrix(NA, nrow = num_leads, ncol = 2)
diff_ci <- matrix(NA, nrow = num_leads, ncol = 2)



treated_BS_ci <- matrix(NA, nrow = num_leads, ncol = 2)
untreated_BS_ci <- matrix(NA, nrow = num_leads, ncol = 2)
diff_BS_ci <- matrix(NA, nrow = num_leads, ncol = 2)

# Confidence level (e.g., 95%)
alpha <- 0.05

# Compute confidence intervals for treated and untreated leads
for (i in 1:num_leads) {
  # Treated
  treated_ci[i, ] <- quantile(bootstrap_replicates[, i], probs = c(alpha / 2, 1 - alpha / 2))
  
  # Untreated
  untreated_ci[i, ] <- quantile(bootstrap_replicates[, i + num_leads], probs = c(alpha / 2, 1 - alpha / 2))
  
  # Treated
  treated_BS_ci[i, ] <- quantile(bootstrap_replicates[, 10 + i], probs = c(alpha / 2, 1 - alpha / 2))
  
  # Untreated
  untreated_BS_ci[i, ] <- quantile(bootstrap_replicates[, 10 + i + num_leads], probs = c(alpha / 2, 1 - alpha / 2))
  
  diff_ci[i, ] <- quantile(bootstrap_replicates[, 20 + i], probs = c(alpha / 2, 1 - alpha / 2))
  diff_BS_ci[i, ] <- quantile(bootstrap_replicates[, 20 + i + num_leads], probs = c(alpha / 2, 1 - alpha / 2))
  
  
}

# Combine results into a data frame for easy visualization
ci_results <- data.frame(
  Lead = 0:(num_leads - 1),
  Treated_Mean = bootstrap_results$t0[1:num_leads],
  Treated_Lower = treated_ci[, 1],
  Treated_Upper = treated_ci[, 2],
  Untreated_Mean = bootstrap_results$t0[(num_leads + 1):(2 * num_leads)],
  Untreated_Lower = untreated_ci[, 1],
  Untreated_Upper = untreated_ci[, 2],
  
  Treated_BS_Mean = bootstrap_results$t0[11:15],
  Treated_BS_Lower = treated_BS_ci[, 1],
  Treated_BS_Upper = treated_BS_ci[, 2],
  Untreated_BS_Mean = bootstrap_results$t0[16:20],
  Untreated_BS_Lower = untreated_BS_ci[, 1],
  Untreated_BS_Upper = untreated_BS_ci[, 2],
  
  Diff_Mean = bootstrap_results$t0[21:25],
  Diff_Lower = diff_ci[, 1],
  Diff_Upper = diff_ci[, 2],
  
  Diff_BS_Mean = bootstrap_results$t0[26:30],
  Diff_BS_Lower = diff_BS_ci[, 1],
  Diff_BS_Upper = diff_BS_ci[, 2]
)

# Print results
print(ci_results)


survivor_plot_boot <- merge(survivor_plot, ci_results, by.x = "time_to_treat", by.y = "Lead", all.x = T)

if (type == "High"){
  starting_BS <- filter(survivor_plot_boot, time_to_treat == -1)$Pct_Burned_HS_Treat
}

survivor_plot_boot[1, c("Treated_Mean", "Treated_Lower", "Treated_Upper", 
                        "Untreated_Mean", "Untreated_Lower", "Untreated_Upper")] <- rep(1, 6)

survivor_plot_boot[1, c("Diff_Mean", "Diff_Lower", "Diff_Upper", 
                        "Diff_BS_Mean", "Diff_BS_Lower", "Diff_BS_Upper")] <- rep(0, 6)


survivor_plot_boot[1, c("Treated_BS_Mean", "Treated_BS_Lower", "Treated_BS_Upper", 
                        "Untreated_BS_Mean", "Untreated_BS_Lower", "Untreated_BS_Upper")] <- rep(starting_BS, 6)

if (type == "High"){
  # Create survivor_BS_plot with no legend
  survival_H_BS_plot_boot <- ggplot(survivor_plot_boot, aes(x = distance)) +
    # Confidence intervals as ribbons
    # Confidence intervals as ribbons
    geom_ribbon(aes(ymin = 100*Treated_BS_Lower, ymax = Treated_BS_Upper*100, fill = "With Treatment (Observed)"), alpha = 0.2) +
    geom_ribbon(aes(ymin = 100*Untreated_BS_Lower, ymax = Untreated_BS_Upper*100, fill = "Without Treatment (Counterfactual)"), alpha = 0.2) +
    # Mean survival lines
    geom_line(aes(y = Pct_Burned_HS_Treat*100, color = "With Treatment (Observed)"), size = 1.2) +
    geom_line(aes(y = Pct_Burned_HS_No_Treat*100, color = "Without Treatment (Counterfactual)"), size = 1.2) +
    
    # Labels and titles
    labs(
      title = "Cumulative High Burn Severity",
      x = "Distance from treatment interaction (km)",
      y = "Predicted high severity burn share (%)",
      color = "Treatment Status",
      fill = "Treatment Status"
    ) +
    
    # Custom colors for lines and ribbons
    scale_color_manual(values = c("With Treatment (Observed)" = "#0072B2", "Without Treatment (Counterfactual)" = "#D55E00")) +
    scale_fill_manual(values = c("With Treatment (Observed)" = "#0072B2", "Without Treatment (Counterfactual)" = "#D55E00")) +
    
    # Theme adjustments
    theme_minimal(base_family = "Helvetica") +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.title = element_text(size = 8, face = "bold"),
      axis.text  = element_text(size = 7),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size  = 6,      # smaller than the y-axis labels
        lineheight = 0.9
      ),
      plot.title   = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                  margin = margin(b = 2)),
      axis.ticks = element_line(linewidth = 0.3),
      plot.margin = margin(3, 3, 3, 3)
    )
}

plot1 <- readRDS(here("output", "event_plotBURNSEV_HIGH.rds"))
plot2 <- readRDS(here("output", "event_plotBURNSEV.rds"))

# Combine the two plots and ensure only one legend
final_layout <- (plot1 | survival_H_BS_plot_boot) / (plot2 | survival_BS_plot_boot) +
  plot_layout(guides = "collect") +  # <- collect shared legend
  plot_annotation(tag_levels = "A")  &
  theme(
    plot.tag = element_text(
      family = "Helvetica",
      face   = "bold",
      size   = 10
    ),
    legend.title = element_text(
      family = "Helvetica",
      face   = "bold",
      size   = 7
    ),
    legend.text = element_text(
      family = "Helvetica",
      size   = 7
    ),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "bottom"
  )

ggsave(here("output", "figures", "EventStudy_Survival_High_Ordinal_BS.pdf"), final_layout, 
       width  = 7.24,          # Science 3-column width
       height = 6,    # ~4.5 in
       units  = "in",
       dpi    = 300)


##################    Benefit-Cost Analysis   ##################

#### Counterfactual Acres Burned - With and Without Treatment

## Option 1 - Get predicted burn probabilities for each grid - multiply it by the acres in a plot 

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

acres_burned_savings <- sum(treated_data$Acres_Burned_No_Treat) - sum(treated_data$Acres_Burned_With_Treat) # Savings only considering 5 leads
acres_burned_savings_pct <- (sum(treated_data$Acres_Burned_No_Treat) - sum(treated_data$Acres_Burned_With_Treat))/sum(treated_data$Acres_Burned_Realized)

acres_burned_savings
acres_burned_savings_pct


acres_burned_MS_HS_savings <- sum(treated_data$Acres_Burned_MS_HS_Without_Treat) - sum(treated_data$Acres_Burned_MS_HS_With_Treat) # Savings only considering 5 leads
acres_burned_MS_HS_savings_pct <- (sum(treated_data$Acres_Burned_MS_HS_Without_Treat) - sum(treated_data$Acres_Burned_MS_HS_With_Treat))/sum(treated_data$Acres_Burned_MS_HS_Realized)

acres_burned_MS_HS_savings
acres_burned_MS_HS_savings_pct

acres_burned_HS_savings <- sum(treated_data$Acres_Burned_HS_Without_Treat) - sum(treated_data$Acres_Burned_HS_With_Treat) # Savings only considering 5 leads
acres_burned_HS_savings_pct <- (sum(treated_data$Acres_Burned_HS_Without_Treat) - sum(treated_data$Acres_Burned_HS_With_Treat))/sum(treated_data$Acres_Burned_HS_Realized)

acres_burned_HS_savings
acres_burned_HS_savings_pct

treated_data <- mutate(treated_data, Structure_Burned_With_Treat = survival_with_treatment*Struc_Count,
                       Structure_Burned_No_Treat = survival_without_treatment*Struc_Count,
                       Structure_Burned_Realized = BURN*Struc_Count, 
                       Structure_Value_Lost_With_Treat = survival_with_treatment*TOT_STRUC_VAL,
                       Structure_Value_Lost_No_Treat = survival_without_treatment*TOT_STRUC_VAL,
                       Structure_Value_Lost_Realized = BURN*TOT_STRUC_VAL)

sum(treated_data$Structure_Burned_With_Treat)
sum(treated_data$Structure_Burned_No_Treat)
sum(treated_data$Structure_Burned_Realized)

struc_savings <- sum(treated_data$Structure_Burned_No_Treat) - sum(treated_data$Structure_Burned_With_Treat) 
struc_savings # Savings only considering 5 leads


(sum(treated_data$Structure_Burned_No_Treat) - sum(treated_data$Structure_Burned_With_Treat))/sum(treated_data$Structure_Burned_Realized)

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
PM25_savings_pct <- (sum(treated_data$PM25_No_Treat, na.rm = T) - sum(treated_data$PM25_With_Treat, na.rm = T))/sum(treated_data$PM25_Realized)

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

# Comment if don't want to use conservative imputed costs

# treated_data <- mutate(treated_data,
#                        Deaths_With_Treat = ifelse(exposure_predicted == T, NA, Acres_Burned_With_Treat*Deaths_PER_ACRE),
#                        Deaths_No_Treat = ifelse(exposure_predicted == T, NA, Acres_Burned_No_Treat*Deaths_PER_ACRE),
#                        Deaths_Realized = ifelse(exposure_predicted == T, NA, Acres_Burned_Realized*Deaths_PER_ACRE),
#                        DeathCost_With_Treat = ifelse(exposure_predicted == T, NA, Acres_Burned_With_Treat*DeathCost_PER_ACRE),
#                        DeathCost_No_Treat = ifelse(exposure_predicted == T, NA, Acres_Burned_No_Treat*DeathCost_PER_ACRE),
#                        DeathCost_Realized = ifelse(exposure_predicted == T, NA, Acres_Burned_Realized*DeathCost_PER_ACRE),
#                        EarningsLoss_With_Treat = ifelse(exposure_predicted == T, NA, Acres_Burned_With_Treat*EarningsLoss_PER_ACRE),
#                        EarningsLoss_No_Treat = ifelse(exposure_predicted == T, NA, Acres_Burned_No_Treat*EarningsLoss_PER_ACRE),
#                        EarningsLoss_Realized = ifelse(exposure_predicted == T, NA, Acres_Burned_Realized*EarningsLoss_PER_ACRE))


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

##### Now lets add additional fires

# Grids_2006_2017 <- read_csv(here("data", "intermediate", "Extra_Fires_v3_24L.csv"))
# 
# extra_treated_data <- Grids_2006_2017 %>% filter(Treated_Dir == 1 & Burn_Right_Away == 0) %>% filter(time_to_treat %in% seq(-1,4,1))
# 
# extra_treated_data <- mutate(extra_treated_data,
#                              Acres_Burned_With_Treat = BURN*Grid_Acres)
# 
# 
# extra_treated_data <- mutate(extra_treated_data, 
#                              PM25_With_Treat = Acres_Burned_With_Treat*FIRE_PM25_PER_ACRE,
#                              CO2_With_Treat = Acres_Burned_With_Treat*FIRE_CO2_PER_ACRE,
#                              Deaths_With_Treat = Acres_Burned_With_Treat*Deaths_PER_ACRE,
#                              DeathCost_With_Treat = Acres_Burned_With_Treat*DeathCost_PER_ACRE,
#                              EarningsLoss_With_Treat = Acres_Burned_With_Treat*EarningsLoss_PER_ACRE)
# 
# 
# acres_burned_extra_savings <- acres_burned_savings_pct*sum(extra_treated_data$Acres_Burned_With_Treat, na.rm = T)
# PM25_extra_savings <- PM25_savings_pct*sum(extra_treated_data$PM25_With_Treat, na.rm = T)
# CO2_extra_savings <- CO2_savings_pct*185*sum(extra_treated_data$CO2_With_Treat, na.rm = T)
# Saved_Deaths_extra <- Saved_Deaths_pct*sum(extra_treated_data$Deaths_With_Treat, na.rm = T)
# DeathCosts_extra_savings <- Death_savings_pct*sum(extra_treated_data$DeathCost_With_Treat, na.rm = T)
# EarningsLoss_extra_savings <- Earnings_savings_pct*sum(extra_treated_data$EarningsLoss_With_Treat, na.rm = T)
# 
# 
# ## Update savings with extra fires
# 
# acres_burned_savings <- acres_burned_savings + acres_burned_extra_savings
# PM25_savings <- PM25_savings + PM25_extra_savings
# carbon_savings <- carbon_savings + CO2_extra_savings
# Saved_Deaths <- Saved_Deaths + Saved_Deaths_extra
# 
# Health_Savings <- Health_Savings + DeathCosts_extra_savings + EarningsLoss_extra_savings

#### Calculate Treatment Costs

# Load Functions 

source(here("code", "functions","tidy_facts.R"))
source(here("code", "functions","calculate_distance.R"))
st_erase = function(x, y) st_difference(x, st_make_valid(st_union(st_combine(y)))) # Taken from here https://r-spatial.github.io/sf/reference/geos_binary_ops.html

########    Load in Datasets    ########

#### FACTS - Forest Service Fuel Treatments

facts <- st_read(here("data", "raw", "FACTS", "S_USA.Activity_HazFuelTrt_PL.shp"))
facts <- tidy_facts(facts, inf_yr = 2023, crs = 5070)
facts$ROW_ID <- 1:nrow(facts)

# facts_filt <- filter(facts, COMPLETE == 1) %>% filter(YEAR >= 2007 & YEAR <= 2013)

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


# #### Option 1 - Calculate lambda by 2006-2013 treatments
# 
# facts_act_grouped <- facts_west %>%
#   filter(YEAR >= 2006 & YEAR <= 2013) %>%
#   group_by(ACTIVITY_UNIT_CN) %>%
#   summarise(geometry = st_union(geometry), 
#             YEAR_COMP = max(YEAR))
# 
# facts_act_grouped$Acres_Treated <- as.numeric(st_area(facts_act_grouped)/4046.86)
# 
# facts_act_grouped <- filter(facts_act_grouped, YEAR_COMP >= 2004)
# 
# 
# #### Intersect FACTS with MTBS to get estimate of lambda
# 
# facts_mtbs_int <- st_intersection(facts_act_grouped, mtbs_filt)
# 
# activities_int <- as.data.frame(facts_mtbs_int) %>%
#   filter(YEAR_COMP <= YEAR_MTBS) %>%
#   filter(YEAR_MTBS - YEAR_COMP <= 10) %>%
#   group_by(ACTIVITY_UNIT_CN) %>%
#   summarise(INTERSECTS = 1)
# 
# facts_act_grouped_df <- merge(as.data.frame(facts_act_grouped), activities_int, by = "ACTIVITY_UNIT_CN",all.x = T)
# facts_act_grouped_df <- dplyr::select(facts_act_grouped_df, -geometry)
# facts_act_grouped_df[is.na(facts_act_grouped_df$INTERSECTS), "INTERSECTS"] <- 0
# 
# treated_data_filt <- filter(treated_data, YEAR %in% seq(2017,2023,1))
# 
# quantiles <- quantile(treated_data_filt$TREAT_SIZE, probs = c(1/3,2/3,1)) 
# 
# q1 <- as.numeric(quantiles[1])
# q2 <- as.numeric(quantiles[2])
# q3 <- as.numeric(quantiles[3])
# 
# lambda_1 <- mean(filter(facts_act_grouped_df, Acres_Treated < q1)$INTERSECTS, na.rm = T)
# lambda_2 <- mean(filter(facts_act_grouped_df, Acres_Treated >= q1 & Acres_Treated <= q2)$INTERSECTS, na.rm = T)
# lambda_3 <- mean(filter(facts_act_grouped_df, Acres_Treated > q2)$INTERSECTS, na.rm = T)
# 
# 
# treated_data_filt <- mutate(treated_data_filt, lambda = case_when(TREAT_SIZE < q1 ~ lambda_1,
#                                                                   TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ lambda_2,
#                                                                   TREAT_SIZE > q2 ~ lambda_3,))
# 
# treated_data_filt <- mutate(treated_data_filt, Size_Bin = case_when(TREAT_SIZE < q1 ~ "Small",
#                                                                     TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ "Medium",
#                                                                     TREAT_SIZE > q2 ~ "Large",))
# 
# treat_ids <- unique(treated_data_filt$TREAT_ID)
# treat_ids <- treat_ids[!is.na(treat_ids)]
# 
# treatment_costs <- treated_data_filt %>%
#   group_by(TREAT_ID) %>%
#   summarise(TREAT_COST = dplyr::first(TREAT_COST), 
#             lambda = dplyr::first(lambda),
#             Size_Bin = dplyr::first(Size_Bin),
#             YEAR_COMPLETE = dplyr::first(YEAR) - dplyr::first(YEARS_SINCE_TREAT) 
#             )
# 
# Benefit_Cost_Ratio_Conditional = (house_val_savings + carbon_savings + Health_Savings)/Treat_Cost


#### Option 2 - Estimate lambda using survival function approach - Kaplan-Meier Estimator

facts_act_grouped <- facts_west_filt %>%
  group_by(ACTIVITY_UNIT_CN) %>%
  summarise(geometry = st_union(geometry), 
            YEAR_COMP = max(YEAR), 
            YEAR_COMP_MIN = min(YEAR))

facts_act_grouped$Acres_Treated <- as.numeric(st_area(facts_act_grouped)/4046.86)

facts_mtbs_int <- st_intersection(facts_act_grouped, mtbs_filt)

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

# Fit Kaplan-Meier survival curve
km_fit <- survfit(Surv(event_time, INTERSECTS) ~ 1, data = facts_act_grouped)

# Plot survival curve
ggsurvplot(km_fit, conf.int = TRUE, ggtheme = theme_minimal())

# Estimate lambda (probability of intersection within 10 years)
lambda_estimate <- 1 - summary(km_fit, times = 10)$surv
print(lambda_estimate)


treated_data_filt <- filter(treated_data, YEAR %in% seq(2017,2023,1))

quantiles <- quantile(treated_data_filt$TREAT_SIZE, probs = c(1/3,2/3,1)) 

# q1 <- as.numeric(quantiles[1])
# q2 <- as.numeric(quantiles[2])
# q3 <- as.numeric(quantiles[3])

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
Avoided_Damages_Baseline <- house_val_savings + carbon_savings + Health_Savings
Benefit_Cost_Ratio_Conditional = Avoided_Damages_Baseline/Treat_Cost



###### Calculate Benefits

Physical_Savings_df <- data.frame(matrix(ncol = 5, nrow = 1))

colnames(Physical_Savings_df) <- c("Acres Burned",  "Structure Loss", "CO2 Emission Savings (t)", "PM2.5 Emissions (t)", "Deaths")

Physical_Savings_df[1, ] <- c(
  prettyNum(round(acres_burned_savings, digits = 0), scientific = F, big.mark = ","),
  prettyNum(round(struc_savings, digits = 0), scientific = F, big.mark = ","),
  prettyNum(round(CO2_savings, digits = 0), scientific = F, big.mark = ","),
  prettyNum(round(PM25_savings, digits = 0), scientific = F, big.mark = ","),
  prettyNum(round(Saved_Deaths, digits = 0), scientific = F, big.mark = ",")
)


Counterfactual_Benefits_df <- data.frame(matrix(ncol = 5, nrow = 1))

colnames(Counterfactual_Benefits_df) <- c("Cost of Treatments",  "Housing Values", "Social Cost of Carbon", "Health Benefits", "CBCR")

Counterfactual_Benefits_df[1, ] <- c(
  paste0("$", prettyNum(round(Treat_Cost, digits = 2), scientific = F, big.mark = ",")),
  paste0("$", prettyNum(round(house_val_savings, digits = 2), scientific = F, big.mark = ",")),
  paste0("$", prettyNum(round(carbon_savings, digits = 2), scientific = F, big.mark = ",")), 
  paste0("$", prettyNum(round(Health_Savings, digits = 2), scientific = F, big.mark = ",")), 
  paste0("$", round(Benefit_Cost_Ratio_Conditional, digits = 2))
)


notes <- "\\parbox[t]{\\textwidth}{\\footnotesize{The following table shows the costs and benefits of fuel treatments}}"

# Convert to a single line string
notes_single_line <- gsub("\n", " ", notes)

stargazer(Counterfactual_Benefits_df, 
          title = ("\\label{tab:CounterfactualAnalysis} Counterfactual Benefits of USFS Fuel Treatments"), 
          summary = F,
          digits = 0,
          type = "latex", 
          out = here("output", "tables", "CounterfactualAnalysis.tex"),
          rownames = F,
          colnames = T, 
          notes = notes_single_line,
          column.sep.width = "-4pt",
          notes.align = "l"
)

stargazer(Physical_Savings_df, 
          title = ("\\label{tab:CounterfactualAnalysis} Counterfactual Benefits of USFS Fuel Treatments"), 
          summary = F,
          digits = 0,
          type = "latex", 
          rownames = F,
          colnames = T, 
          notes = notes_single_line,
          column.sep.width = "-4pt",
          notes.align = "l"
)


###### Calculate Size Specific Benefits

## Cost per size bin

Treat_Size_Benefits <- treated_data_filt %>%
  group_by(Size_Bin) %>%
  summarise(Acres_Burned_Savings = sum(Acres_Burned_No_Treat - Acres_Burned_With_Treat),
            Struc_Savings = sum(Structure_Burned_No_Treat - Structure_Burned_With_Treat),
            House_Savings = sum(Structure_Value_Lost_No_Treat - Structure_Value_Lost_With_Treat, na.rm = T),
            Carbon_Savings = (sum(CO2_No_Treat - CO2_With_Treat))*185,
            Death_Savings = sum(DeathCost_No_Treat - DeathCost_With_Treat, na.rm = T),
            Earning_Savings = sum(EarningsLoss_No_Treat - EarningsLoss_With_Treat, na.rm = T),
            N_Treats = n_distinct(TREAT_ID),
            lambda = dplyr::first(lambda))

Treat_Size_Benefits <- mutate(Treat_Size_Benefits, Total_Benefit = House_Savings + Carbon_Savings + Death_Savings)
Treat_Size_Benefits <- mutate(Treat_Size_Benefits, Benefit_Per_Treat = Total_Benefit/N_Treats)

#### Benefit per Treatment

sum(Treat_Size_Benefits$Total_Benefit)/sum(Treat_Size_Benefits$N_Treats)

#### Benefit Cost Ratios for different sized treatments

mu_1 <- filter(Treat_Size_Benefits, Size_Bin == "Small")$Benefit_Per_Treat
mu_2 <- filter(Treat_Size_Benefits, Size_Bin == "Medium")$Benefit_Per_Treat
mu_3 <- filter(Treat_Size_Benefits, Size_Bin == "Large")$Benefit_Per_Treat

#### Average Benefit Per Acres

Treat_Sizes <- treated_data_filt %>%
  group_by(TREAT_ID) %>%
  summarise(Size_Bin = dplyr::first(Size_Bin), 
            TREAT_SIZE = dplyr::first(TREAT_SIZE),
  )

Size_Bin_Acres_Treated <- Treat_Sizes %>%
  group_by(Size_Bin) %>%
  summarise(Total_Acres_Treated = sum(TREAT_SIZE))

Treat_Size_Benefits <- merge(Treat_Size_Benefits, Size_Bin_Acres_Treated, by = "Size_Bin", all.x = T)

Treat_Size_Benefits <- mutate(Treat_Size_Benefits, Benefit_Per_Acre = Total_Benefit/Total_Acres_Treated)

mu_perA_1 <- filter(Treat_Size_Benefits, Size_Bin == "Small")$Benefit_Per_Acre
mu_perA_2 <- filter(Treat_Size_Benefits, Size_Bin == "Medium")$Benefit_Per_Acre
mu_perA_3 <- filter(Treat_Size_Benefits, Size_Bin == "Large")$Benefit_Per_Acre



#### Costs  based on bin size

Cost_Size_Bin <- treatment_costs %>%
  group_by(Size_Bin) %>%
  summarise(Cost = sum(TREAT_COST_REP, na.rm = T))

C_1 <- filter(Cost_Size_Bin, Size_Bin == "Small")$Cost
C_2 <- filter(Cost_Size_Bin, Size_Bin == "Medium")$Cost
C_3 <- filter(Cost_Size_Bin, Size_Bin == "Large")$Cost

B_1 <- filter(Treat_Size_Benefits, Size_Bin == "Small")$Total_Benefit
B_2 <- filter(Treat_Size_Benefits, Size_Bin == "Medium")$Total_Benefit
B_3 <- filter(Treat_Size_Benefits, Size_Bin == "Large")$Total_Benefit

#### Conditional B-C Ratios & Conditional ROI

BC_1 <- B_1/C_1
BC_2 <- B_2/C_2
BC_3 <- B_3/C_3

ROI_1 <- BC_1^(1/10)
ROI_2 <- BC_2^(1/10)
ROI_3 <- BC_3^(1/10)

c(BC_1, BC_2, BC_3)
c(ROI_1, ROI_2, ROI_3)

#### Calculate the Ex-ante B/C Ratio for treatments completed 2007-2023

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


facts_act_grouped_06_23 <- mutate(facts_act_grouped_06_23, lambda = case_when(TREAT_SIZE < q1 ~ lambda_1,
                                                                              TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ lambda_2,
                                                                              TREAT_SIZE > q2 ~ lambda_3),
                                  mu = case_when(TREAT_SIZE < q1 ~ mu_1,
                                                 TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ mu_2,
                                                 TREAT_SIZE > q2 ~ mu_3))

facts_act_grouped_06_23 <- mutate(facts_act_grouped_06_23, Size_Bin = case_when(TREAT_SIZE < q1 ~ "Small",
                                                                                TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ "Medium",
                                                                                TREAT_SIZE > q2 ~ "Large"))


facts_act_grouped_06_23 <- mutate(facts_act_grouped_06_23, Benefits = mu*lambda, 
                                  Benefit_Cost = Benefits/TOT_COST)

# facts_act_grouped_06_23 <- filter(facts_act_grouped_06_23, ACCOMPLISHED_ACRES > 0)


## Identify observations with high cost/acre

mean_cost <- mean(facts_west$COST_PER_UOM, na.rm = TRUE)
sd_cost <- sd(facts_west$COST_PER_UOM, na.rm = TRUE)
threshold <- mean_cost + 10*sd_cost


facts_act_grouped_06_23_filt <- filter(facts_act_grouped_06_23, COST_PER_UOM <= threshold)

ex <- filter(facts_act_grouped_06_23, ACTIVITY_UNIT_CN == "6150108010602") # one treatment that cost 1 billion USD

BCR_Overall <- sum(facts_act_grouped_06_23_filt$Benefits)/sum(facts_act_grouped_06_23_filt$TOT_COST)

BCR_Overall

mean(facts_act_grouped_06_23_filt$Benefit_Cost)
median(facts_act_grouped_06_23_filt$Benefit_Cost)

facts_act_grouped_06_23_filt_1 <- filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q_05)

BCR_Overall <- sum(facts_act_grouped_06_23_filt_1$Benefits)/sum(facts_act_grouped_06_23_filt_1$TOT_COST)
mean(facts_act_grouped_06_23_filt_1$Benefit_Cost)
Median_BC <- median(facts_act_grouped_06_23_filt_1$Benefit_Cost)


#### Ex ante B-C Ratio by Treatment Size Bins: "Small" (75-600), "Medium" (600-2400), "Large" (>2400)

q_05 <- quantile(treated_data_filt$TREAT_SIZE, probs = .05) 

C_1 <- mean(filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q_05 & TREAT_SIZE <= q1)$TOT_COST, na.rm = T)
C_2 <- mean(filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q1 & TREAT_SIZE <= q2)$TOT_COST, na.rm = T)
C_3 <- mean(filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q2)$TOT_COST, na.rm = T)

BC_Large <- sum(filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q2)$Benefits)/sum(filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q2)$TOT_COST)
BC_Medium <- sum(filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q1 & TREAT_SIZE <= q2)$Benefits)/sum(filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q1 & TREAT_SIZE <= q2)$TOT_COST)
BC_Small <- sum(filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q_05 & TREAT_SIZE <= q1)$Benefits)/sum(filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q_05 & TREAT_SIZE <= q1)$TOT_COST)


c(BC_Large, BC_Medium, BC_Small)


## Create Table of Robustness of BCR Estimates by different specifications

N_Treats_Baseline <- sum(Treat_Size_Benefits$N_Treats)

Robustness <- data.frame(matrix(ncol = 3, nrow = 5))

colnames(Robustness) <- c("Ex-Ante BCR", "No. Treatments", "Sample Period")
rownames(Robustness) <- c("Baseline",  "Alternative Cost Estimate", "No Smoke Exposure Imputation", 
                          "2017-2020 Fires No Imputation", "2017-2020 Fires Imputation")

Robustness[1, ] <- c(
  round(BCR_Overall, digits = 2), 
  N_Treats_Baseline, 
  "2017-2023"
)

###### Alternative Benefit-Cost Ratio's with estimated costs

AVG_COST_FOOTPRINT_ACRE <- 400.9842 # taken from "DescriptiveStats.R"

facts_act_grouped_06_23 <- mutate(facts_act_grouped_06_23, TOT_COST_ALT = TREAT_SIZE * AVG_COST_FOOTPRINT_ACRE,
                                  Benefit_Cost_Alt = Benefits/TOT_COST_ALT)

facts_act_grouped_06_23_filt_2 <- filter(facts_act_grouped_06_23, TREAT_SIZE > q_05)

BCR_Overall_Alt <- sum(facts_act_grouped_06_23_filt_2$Benefits)/sum(facts_act_grouped_06_23_filt_2$TOT_COST_ALT)
BCR_Overall_Alt
BC_Large_Alt <- sum(filter(facts_act_grouped_06_23_filt_2, TREAT_SIZE > q2)$Benefits)/sum(filter(facts_act_grouped_06_23_filt_2, TREAT_SIZE > q2)$TOT_COST_ALT)
BC_Medium_Alt <- sum(filter(facts_act_grouped_06_23_filt_2, TREAT_SIZE > q1 & TREAT_SIZE <= q2)$Benefits)/sum(filter(facts_act_grouped_06_23_filt_2, TREAT_SIZE > q1 & TREAT_SIZE <= q2)$TOT_COST_ALT)
BC_Small_Alt <- sum(filter(facts_act_grouped_06_23_filt_2, TREAT_SIZE > q_05 & TREAT_SIZE <= q1)$Benefits)/sum(filter(facts_act_grouped_06_23_filt_2, TREAT_SIZE > q_05 & TREAT_SIZE <= q1)$TOT_COST_ALT)


c(BC_Large_Alt, BC_Medium_Alt, BC_Small_Alt)

Robustness[2, ] <- c(
  round(BCR_Overall_Alt, digits = 2),
  N_Treats_Baseline,
  "2017-2023"
)


###### Calculate BCRs with discounting

# Define discount rates and years
discount_rates <- c(0.03, 0.05, 0.08)
years <- 0:10

# Create a grid of years and discount rates for expansion
discount_df <- expand.grid(t = years, r = discount_rates)

# Add yearly discounted benefits to each treatment row
facts_discounted <- facts_act_grouped_06_23_filt %>%
  dplyr::select(ACTIVITY_UNIT_CN, TREAT_SIZE, lambda, mu, TOT_COST) %>%
  crossing(discount_df) %>%
  mutate(Benefit_t = mu * (lambda / 10) * (1 - lambda / 10)^t * (1 + r)^(-t))


# Sum discounted benefits across time (NPV)
facts_npv <- facts_discounted %>%
  group_by(ACTIVITY_UNIT_CN, TREAT_SIZE, TOT_COST, r) %>%
  summarise(NPV_Benefit = sum(Benefit_t), .groups = "drop")

# Define quantile cutoffs (already assumed defined: q_05, q1, q2)

# Classify treatment sizes
facts_npv <- facts_npv %>%
  mutate(Size = case_when(
    TREAT_SIZE > q2 ~ "Large",
    TREAT_SIZE > q1 & TREAT_SIZE <= q2 ~ "Medium",
    TREAT_SIZE > q_05 & TREAT_SIZE <= q1 ~ "Small",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Size))

# Compute B-C ratios by size and discount rate
bc_by_size <- facts_npv %>%
  group_by(r, Size) %>%
  summarise(
    Total_Benefits = sum(NPV_Benefit),
    Total_Costs = sum(TOT_COST),
    BC_Ratio = Total_Benefits / Total_Costs,
    .groups = "drop"
  ) %>%
  dplyr::select(r, Size, BC_Ratio) %>%
  pivot_wider(names_from = Size, values_from = BC_Ratio) %>%
  arrange(r)

# Compute overall B-C ratio (all sizes)
bc_total <- facts_npv %>%
  group_by(r) %>%
  summarise(
    Total_Benefits = sum(NPV_Benefit),
    Total_Costs = sum(TOT_COST),
    Total = Total_Benefits / Total_Costs,
    .groups = "drop"
  )

# Join and reorder columns
bc_table <- left_join(bc_by_size, bc_total, by = "r") %>%
  dplyr::select(r, Small, Medium, Large, Total)

# Optional: round for display
bc_table <- bc_table %>%
  mutate(across(where(is.numeric), round, 2))

bc_table <- rbind(round(c(0, BC_Small, BC_Medium, BC_Large, BCR_Overall), digits = 2), bc_table)

# Optional: format for display
bc_table <- mutate(bc_table, across(where(is.numeric), round, 2))

bc_table

colnames(bc_table) <- c("Discount Rate", "Small", "Medium", "Large", "Total")

notes <- "\\parbox[t]{\\textwidth}{\\footnotesize{Estimated benefit–cost ratios for small (75–600 acres), medium (600–2400 acres), and large (>2400 acres) treatments, as well as the overall benefit–cost ratio, reported by discount rate.}}"

# Convert to a single line string
notes_single_line <- gsub("\n", " ", notes)


stargazer(bc_table, 
          title = ("\\label{tab:BenefitCostRatios} Estimated Benefit-Cost Ratios by Discount Rate"), 
          summary = F,
          digits = 2,
          type = "latex", 
          out = here("output", "tables", "BenefitCostRatios.tex"),
          rownames = F,
          colnames = T, 
          notes = notes_single_line,
          column.sep.width = "-4pt",
          notes.align = "l"
)


#### Calculate Ex Ante Annualized Rate of Return

## ROI = [(mu * lambda)/C]^[1/T], mu is average yearly benefits, lambda is probability of intersecting a fire in 10 years, C is cost of a fuel treatment

ROI_Small <- BC_Small^(1/10)
ROI_Medium <- BC_Medium^(1/10)
ROI_Large <- BC_Large^(1/10)

c(ROI_Small, ROI_Medium, ROI_Large)


########   Distribution of lambda and mu by treatment size

#### Survival plots for each treatment size

# Add a treatment size category
facts_act_grouped <- facts_act_grouped %>%
  mutate(Treatment_Size = case_when(
    Acres_Treated > q_05 & Acres_Treated < q1 ~ "Small",
    Acres_Treated >= q1 & Acres_Treated <= q2 ~ "Medium",
    Acres_Treated > q2 ~ "Large"
  ))

# Define professional color palette
nature_palette <- c("#0072B2", "#D55E00", "#33a02c") # Blue, Green, Red

# Fit Kaplan-Meier survival curves
km_fit <- survfit(Surv(event_time, INTERSECTS) ~ Treatment_Size, data = facts_act_grouped)

# Extract survival estimates into dataframe
km_df <- data.frame(
  time = km_fit$time,
  surv = km_fit$surv,
  strata = rep(names(km_fit$strata), km_fit$strata),
  lower = km_fit$lower,  # Confidence intervals (optional)
  upper = km_fit$upper
)

# Ensure all curves start at (0,1)
start_df <- km_df %>%
  group_by(strata) %>%
  summarise(time = 0, surv = 1, .groups = "drop") 

km_df <- bind_rows(start_df, km_df) %>% arrange(strata, time)

# Add vertical drop segments
km_df <- km_df %>%
  group_by(strata) %>%
  mutate(
    time_lag = lag(time, default = 0),  # Previous time point
    surv_lag = lag(surv, default = 1)   # Previous survival value
  )

# Remove the prefix from strata labels
km_df$strata <- gsub("Treatment_Size=", "", km_df$strata)

# Ensure correct ordering of factor levels
km_df$strata <- factor(km_df$strata, levels = c("Small", "Medium", "Large"))

# Plot Kaplan-Meier with vertical drop lines
lambda_plot <- ggplot(km_df, aes(x = time, y = surv, color = strata)) +
  geom_step(size = 1.2, alpha = 2) +  # Main survival curve
  geom_segment(aes(x = time, xend = time, y = surv_lag, yend = surv), linetype = "dashed", size = 0.8) +  # Vertical drops
  scale_color_manual(values = nature_palette, labels = c("Small (75-600)", "Medium (600-2400)", "Large (> 2400)")) +
  labs(
    x = "Time since treatment (Years)",
    y = "Probability of no fire interaction",
    color = NULL  # Removes default "strata" title from legend
  ) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text  = element_text(size = 7),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size  = 6,      # smaller than the y-axis labels
      lineheight = 0.9
    ),
    plot.title   = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                margin = margin(b = 2)),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(3, 3, 3, 3)
  ) +
  coord_cartesian(clip = "off")  # Prevents cropping of axis lines

final_plot <- lambda_plot +
  plot_layout(guides = "collect") &   # <- collect shared legend
  theme(legend.position = "bottom")  # <- move legend to bottom


ggsave(here("output", "figures", "TreatmentInteractionSurvival.pdf"), final_plot, 
       width  = 7.24,          # Science 3-column width
       height = 4.5,    # ~4.22 in
       units  = "in",
       dpi    = 300)



#### Survival Plots for each treatment size

q1 <- 600
q2 <- 2400

type <- "Med-High"

large_plot <- calculate_size_survival("Large")
medium_plot <- calculate_size_survival("Medium")
small_plot <- calculate_size_survival("Small")

large_survival_plot <- large_plot$survivor_plot
medium_survival_plot <- medium_plot$survivor_plot
small_survival_plot <- small_plot$survivor_plot

large_survival_BS_plot <- large_plot$survivor_BS_plot
medium_survival_BS_plot <- medium_plot$survivor_BS_plot
small_survival_BS_plot <- small_plot$survivor_BS_plot

# Combine the two plots and ensure only one legend
final_layout <- (small_survival_plot | medium_survival_plot | large_survival_plot) /
  (small_survival_BS_plot | medium_survival_BS_plot | large_survival_BS_plot) +
  plot_layout(guides = "collect") +  # <- collect shared legend
  plot_annotation(tag_levels = "A")  &
  theme(
    plot.tag = element_text(
      family = "Helvetica",
      face   = "bold",
      size   = 10
    ),
    legend.title = element_text(
      family = "Helvetica",
      face   = "bold",
      size   = 7
    ),
    legend.text = element_text(
      family = "Helvetica",
      size   = 7
    ),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "bottom",
    plot.margin = margin(3, 10, 3, 3)  # <-- add right margin (try 8–14)
    
  )

ggsave(here("output", "figures", "SurvivalSizeHeterogeneity.pdf"), final_layout,
       width  = 7.24,          # Science 3-column width
       height = 8,    # ~4.22 in
       units  = "in",
       dpi    = 300)



row_title1 <- wrap_elements(
  full = grid::textGrob(
    "Unconditional Burn Probability", 
    gp = gpar(fontsize = 20, fontface = "bold")
  ))

row_title2 <- wrap_elements(
  full = grid::textGrob(
    "Unconditional Burn Severity", 
    gp = gpar(fontsize = 20, fontface = "bold")
  ))

# Combine plots into rows
burn_prob_row <- small_survival_plot | medium_survival_plot | large_survival_plot
burn_sev_row  <- small_survival_BS_plot | medium_survival_BS_plot | large_survival_BS_plot

# Final layout with row titles
final_layout <- (
  burn_prob_row / burn_sev_row
) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a",
                  title = "Unconditional Burn Probability & Burn Severity by Treatment Size Category",
                  theme = theme(
                    plot.title = element_text(size = 22, face = "bold"),
                    plot.subtitle = element_text(size = 12)
                  )
  ) &
  theme(legend.position = "bottom")

# Print the final plot
final_layout


ggsave(here("output", "figures", "SurvivalSizeHeterogeneity.pdf"), final_layout, width = 16, height = 12, units = "in")




###### Drop Imputation Smoke costs in BCR calculation (2017-2023)  ######

treated_data <- mutate(treated_data,
                       Deaths_With_Treat = ifelse(exposure_predicted == T, NA, Acres_Burned_With_Treat*Deaths_PER_ACRE),
                       Deaths_No_Treat = ifelse(exposure_predicted == T, NA, Acres_Burned_No_Treat*Deaths_PER_ACRE),
                       Deaths_Realized = ifelse(exposure_predicted == T, NA, Acres_Burned_Realized*Deaths_PER_ACRE),
                       DeathCost_With_Treat = ifelse(exposure_predicted == T, NA, Acres_Burned_With_Treat*DeathCost_PER_ACRE),
                       DeathCost_No_Treat = ifelse(exposure_predicted == T, NA, Acres_Burned_No_Treat*DeathCost_PER_ACRE),
                       DeathCost_Realized = ifelse(exposure_predicted == T, NA, Acres_Burned_Realized*DeathCost_PER_ACRE),
                       EarningsLoss_With_Treat = ifelse(exposure_predicted == T, NA, Acres_Burned_With_Treat*EarningsLoss_PER_ACRE),
                       EarningsLoss_No_Treat = ifelse(exposure_predicted == T, NA, Acres_Burned_No_Treat*EarningsLoss_PER_ACRE),
                       EarningsLoss_Realized = ifelse(exposure_predicted == T, NA, Acres_Burned_Realized*EarningsLoss_PER_ACRE))


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

treated_data_filt <- filter(treated_data, YEAR %in% seq(2017,2023,1))

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
Avoided_Damages_Baseline <- house_val_savings + carbon_savings + Health_Savings
Benefit_Cost_Ratio_Conditional = Avoided_Damages_Baseline/Treat_Cost

Treat_Size_Benefits <- treated_data_filt %>%
  group_by(Size_Bin) %>%
  summarise(Acres_Burned_Savings = sum(Acres_Burned_No_Treat - Acres_Burned_With_Treat),
            Struc_Savings = sum(Structure_Burned_No_Treat - Structure_Burned_With_Treat),
            House_Savings = sum(Structure_Value_Lost_No_Treat - Structure_Value_Lost_With_Treat, na.rm = T),
            Carbon_Savings = (sum(CO2_No_Treat - CO2_With_Treat))*185,
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

#### Calculate the Ex-ante B/C Ratio for treatments completed 2007-2023

facts_act_grouped_06_23 <- mutate(facts_act_grouped_06_23, lambda = case_when(TREAT_SIZE < q1 ~ lambda_1,
                                                                              TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ lambda_2,
                                                                              TREAT_SIZE > q2 ~ lambda_3),
                                  mu = case_when(TREAT_SIZE < q1 ~ mu_1,
                                                 TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ mu_2,
                                                 TREAT_SIZE > q2 ~ mu_3))

facts_act_grouped_06_23 <- mutate(facts_act_grouped_06_23, Size_Bin = case_when(TREAT_SIZE < q1 ~ "Small",
                                                                                TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ "Medium",
                                                                                TREAT_SIZE > q2 ~ "Large"))


facts_act_grouped_06_23 <- mutate(facts_act_grouped_06_23, Benefits = mu*lambda, 
                                  Benefit_Cost = Benefits/TOT_COST)


## Identify observations with high cost/acre

mean_cost <- mean(facts_west$COST_PER_UOM, na.rm = TRUE)
sd_cost <- sd(facts_west$COST_PER_UOM, na.rm = TRUE)
threshold <- mean_cost + 10*sd_cost


facts_act_grouped_06_23_filt <- filter(facts_act_grouped_06_23, COST_PER_UOM <= threshold)
facts_act_grouped_06_23_filt_1 <- filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q_05)

BCR_Overall_noimp_2017_2023 <- sum(facts_act_grouped_06_23_filt_1$Benefits)/sum(facts_act_grouped_06_23_filt_1$TOT_COST)

#### Ex ante B-C Ratio by Treatment Size Bins: "Small" (75-600), "Medium" (600-2400), "Large" (>2400)

N_Treats_noimp_2017_2023 <- sum(Treat_Size_Benefits$N_Treats)

Robustness[3, ] <- c(
  round(BCR_Overall_noimp_2017_2023, digits = 2),
  N_Treats_noimp_2017_2023,
  "2017-2023"
)


###### No Imputation Smoke costs for BCR calculation (2017-2020)  ######

treated_data_noimp <- filter(treated_data, exposure_predicted == F) # filter sample to only non-imputed smoke benefits

treated_data_noimp <- mutate(treated_data_noimp, Acres_Burned_With_Treat = survival_with_treatment*Grid_Acres,
                             Acres_Burned_No_Treat = survival_without_treatment*Grid_Acres,
                             Acres_Burned_Realized = BURN*Grid_Acres,
                             Acres_Burned_MS_HS_With_Treat = survival_with_treatment*Grid_Acres*pred_with_treatment_MED_HIGH_BS, 
                             Acres_Burned_MS_HS_Without_Treat = survival_without_treatment*Grid_Acres*pred_without_treatment_MED_HIGH_BS,
                             Acres_Burned_MS_HS_Realized = BURN*Grid_Acres*BURN_SEV_MED_HIGH_PCT,
                             Acres_Burned_HS_With_Treat = survival_with_treatment*Grid_Acres*pred_with_treatment_HIGH_BS, 
                             Acres_Burned_HS_Without_Treat = survival_without_treatment*Grid_Acres*pred_without_treatment_HIGH_BS,
                             Acres_Burned_HS_Realized = BURN*Grid_Acres*BURN_SEV_HIGH_PCT
)

treated_data_noimp <- mutate(treated_data_noimp, Structure_Burned_With_Treat = survival_with_treatment*Struc_Count,
                             Structure_Burned_No_Treat = survival_without_treatment*Struc_Count,
                             Structure_Burned_Realized = BURN*Struc_Count, 
                             Structure_Value_Lost_With_Treat = survival_with_treatment*TOT_STRUC_VAL,
                             Structure_Value_Lost_No_Treat = survival_without_treatment*TOT_STRUC_VAL,
                             Structure_Value_Lost_Realized = BURN*TOT_STRUC_VAL)

carbon_savings <- (sum(treated_data_noimp$CO2_No_Treat, na.rm = T) - sum(treated_data_noimp$CO2_With_Treat, na.rm = T))*185 # CO2 * Carbon price
carbon_savings

#### Structure Value Savings

house_val_savings <- sum(treated_data_noimp$Structure_Value_Lost_No_Treat, na.rm = T) - sum(treated_data_noimp$Structure_Value_Lost_With_Treat, na.rm = T) 
house_val_savings

#### Smoke Emissions

treated_data_noimp <- mutate(treated_data_noimp, PM25_With_Treat = Acres_Burned_With_Treat*FIRE_PM25_PER_ACRE,
                             PM25_No_Treat = Acres_Burned_No_Treat*FIRE_PM25_PER_ACRE,
                             PM25_Realized = Acres_Burned_Realized*FIRE_PM25_PER_ACRE,
                             CO2_With_Treat = Acres_Burned_With_Treat*FIRE_CO2_PER_ACRE,
                             CO2_No_Treat = Acres_Burned_No_Treat*FIRE_CO2_PER_ACRE,
                             CO2_Realized = Acres_Burned_Realized*FIRE_CO2_PER_ACRE)

#### Health Costs

treated_data_noimp <- mutate(treated_data_noimp, 
                             Deaths_With_Treat = Acres_Burned_With_Treat*Deaths_PER_ACRE,
                             Deaths_No_Treat = Acres_Burned_No_Treat*Deaths_PER_ACRE,
                             Deaths_Realized = Acres_Burned_Realized*Deaths_PER_ACRE,
                             DeathCost_With_Treat = Acres_Burned_With_Treat*DeathCost_PER_ACRE,
                             DeathCost_No_Treat = Acres_Burned_No_Treat*DeathCost_PER_ACRE,
                             DeathCost_Realized = Acres_Burned_Realized*DeathCost_PER_ACRE,
                             EarningsLoss_With_Treat = Acres_Burned_With_Treat*EarningsLoss_PER_ACRE,
                             EarningsLoss_No_Treat = Acres_Burned_No_Treat*EarningsLoss_PER_ACRE,
                             EarningsLoss_Realized = Acres_Burned_Realized*EarningsLoss_PER_ACRE)

Death_savings <- (sum(treated_data_noimp$DeathCost_No_Treat, na.rm = T) - sum(treated_data_noimp$DeathCost_With_Treat, na.rm = T)) # Deaths*VSL
Death_savings_pct <- (sum(treated_data_noimp$DeathCost_No_Treat, na.rm = T) - sum(treated_data_noimp$DeathCost_With_Treat, na.rm = T))/sum(treated_data$DeathCost_Realized, na.rm = T)

Death_savings

Earnings_savings <- (sum(treated_data_noimp$EarningsLoss_No_Treat, na.rm = T) - sum(treated_data_noimp$EarningsLoss_With_Treat, na.rm = T)) # Earnings Estimate
Earnings_savings_pct <- (sum(treated_data_noimp$EarningsLoss_No_Treat, na.rm = T) - sum(treated_data_noimp$EarningsLoss_With_Treat, na.rm = T))/sum(treated_data_noimp$EarningsLoss_Realized, na.rm = T)

Earnings_savings

Health_Savings_true <- Death_savings + Earnings_savings

Avoided_Damages_Noimp <- house_val_savings + carbon_savings + Health_Savings_true

treated_data_noimp <- mutate(treated_data_noimp, lambda = case_when(TREAT_SIZE < q1 ~ lambda_1,
                                                                    TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ lambda_2,
                                                                    TREAT_SIZE > q2 ~ lambda_3,))

treated_data_noimp <- mutate(treated_data_noimp, Size_Bin = case_when(TREAT_SIZE < q1 ~ "Small",
                                                                      TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ "Medium",
                                                                      TREAT_SIZE > q2 ~ "Large",))

treat_ids <- unique(treated_data_noimp$TREAT_ID)
treat_ids <- treat_ids[!is.na(treat_ids)]

treatment_costs <- treated_data_noimp %>%
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
Benefit_Cost_Ratio_Conditional_NoImp = Avoided_Damages_Noimp/Treat_Cost

Treat_Size_Benefits_noimp <- treated_data_noimp %>%
  group_by(Size_Bin) %>%
  summarise(Acres_Burned_Savings = sum(Acres_Burned_No_Treat - Acres_Burned_With_Treat),
            Struc_Savings = sum(Structure_Burned_No_Treat - Structure_Burned_With_Treat),
            House_Savings = sum(Structure_Value_Lost_No_Treat - Structure_Value_Lost_With_Treat, na.rm = T),
            Carbon_Savings = (sum(CO2_No_Treat - CO2_With_Treat))*185,
            Death_Savings = sum(DeathCost_No_Treat - DeathCost_With_Treat, na.rm = T),
            Earning_Savings = sum(EarningsLoss_No_Treat - EarningsLoss_With_Treat, na.rm = T),
            N_Treats = n_distinct(TREAT_ID),
            lambda = dplyr::first(lambda))

Treat_Size_Benefits_noimp <- mutate(Treat_Size_Benefits_noimp, Total_Benefit = House_Savings + Carbon_Savings + Death_Savings)
Treat_Size_Benefits_noimp <- mutate(Treat_Size_Benefits_noimp, Benefit_Per_Treat = Total_Benefit/N_Treats)

#### Benefit per Treatment

sum(Treat_Size_Benefits_noimp$Total_Benefit)/sum(Treat_Size_Benefits_noimp$N_Treats)

#### Benefit Cost Ratios for different sized treatments

mu_1_noimp <- filter(Treat_Size_Benefits_noimp, Size_Bin == "Small")$Benefit_Per_Treat
mu_2_noimp <- filter(Treat_Size_Benefits_noimp, Size_Bin == "Medium")$Benefit_Per_Treat
mu_3_noimp <- filter(Treat_Size_Benefits_noimp, Size_Bin == "Large")$Benefit_Per_Treat


#### Calculate the Ex-ante B/C Ratio for treatments completed 2007-2023

facts_act_grouped_06_23 <- mutate(facts_act_grouped_06_23, lambda = case_when(TREAT_SIZE < q1 ~ lambda_1,
                                                                              TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ lambda_2,
                                                                              TREAT_SIZE > q2 ~ lambda_3),
                                  mu = case_when(TREAT_SIZE < q1 ~ mu_1_noimp,
                                                 TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ mu_2_noimp,
                                                 TREAT_SIZE > q2 ~ mu_3_noimp))

facts_act_grouped_06_23 <- mutate(facts_act_grouped_06_23, Size_Bin = case_when(TREAT_SIZE < q1 ~ "Small",
                                                                                TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ "Medium",
                                                                                TREAT_SIZE > q2 ~ "Large"))


facts_act_grouped_06_23 <- mutate(facts_act_grouped_06_23, Benefits = mu*lambda, 
                                  Benefit_Cost = Benefits/TOT_COST)

# facts_act_grouped_06_23 <- filter(facts_act_grouped_06_23, ACCOMPLISHED_ACRES > 0)


## Identify observations with high cost/acre

mean_cost <- mean(facts_west$COST_PER_UOM, na.rm = TRUE)
sd_cost <- sd(facts_west$COST_PER_UOM, na.rm = TRUE)
threshold <- mean_cost + 10*sd_cost


facts_act_grouped_06_23_filt <- filter(facts_act_grouped_06_23, COST_PER_UOM <= threshold)
facts_act_grouped_06_23_filt_1 <- filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q_05)

BCR_Overall_noimp <- sum(facts_act_grouped_06_23_filt_1$Benefits)/sum(facts_act_grouped_06_23_filt_1$TOT_COST)

#### Ex ante B-C Ratio by Treatment Size Bins: "Small" (75-600), "Medium" (600-2400), "Large" (>2400)

N_Treats_noimp <- sum(Treat_Size_Benefits_noimp$N_Treats)

Robustness[4, ] <- c(
  round(BCR_Overall_noimp, digits = 2),
  N_Treats_noimp,
  "2017-2020"
)

###### Imputated Smoke costs (2017-2020 fires) for BCR calculation ######

treated_data_imp <- filter(treated_data, exposure_predicted == F)

treated_data_imp <- mutate(treated_data_imp, Acres_Burned_With_Treat = survival_with_treatment*Grid_Acres,
                           Acres_Burned_No_Treat = survival_without_treatment*Grid_Acres,
                           Acres_Burned_Realized = BURN*Grid_Acres,
                           Acres_Burned_MS_HS_With_Treat = survival_with_treatment*Grid_Acres*pred_with_treatment_MED_HIGH_BS, 
                           Acres_Burned_MS_HS_Without_Treat = survival_without_treatment*Grid_Acres*pred_without_treatment_MED_HIGH_BS,
                           Acres_Burned_MS_HS_Realized = BURN*Grid_Acres*BURN_SEV_MED_HIGH_PCT,
                           Acres_Burned_HS_With_Treat = survival_with_treatment*Grid_Acres*pred_with_treatment_HIGH_BS, 
                           Acres_Burned_HS_Without_Treat = survival_without_treatment*Grid_Acres*pred_without_treatment_HIGH_BS,
                           Acres_Burned_HS_Realized = BURN*Grid_Acres*BURN_SEV_HIGH_PCT
)

treated_data_imp <- mutate(treated_data_imp, Structure_Burned_With_Treat = survival_with_treatment*Struc_Count,
                           Structure_Burned_No_Treat = survival_without_treatment*Struc_Count,
                           Structure_Burned_Realized = BURN*Struc_Count, 
                           Structure_Value_Lost_With_Treat = survival_with_treatment*TOT_STRUC_VAL,
                           Structure_Value_Lost_No_Treat = survival_without_treatment*TOT_STRUC_VAL,
                           Structure_Value_Lost_Realized = BURN*TOT_STRUC_VAL)

carbon_savings <- (sum(treated_data_imp$CO2_No_Treat, na.rm = T) - sum(treated_data_imp$CO2_With_Treat, na.rm = T))*185 # CO2 * Carbon price
carbon_savings

#### Structure Value Savings

house_val_savings <- sum(treated_data_imp$Structure_Value_Lost_No_Treat, na.rm = T) - sum(treated_data_imp$Structure_Value_Lost_With_Treat, na.rm = T) 
house_val_savings

#### Smoke Emissions

treated_data_imp <- mutate(treated_data_imp, PM25_With_Treat = Acres_Burned_With_Treat*FIRE_PM25_PER_ACRE,
                           PM25_No_Treat = Acres_Burned_No_Treat*FIRE_PM25_PER_ACRE,
                           PM25_Realized = Acres_Burned_Realized*FIRE_PM25_PER_ACRE,
                           CO2_With_Treat = Acres_Burned_With_Treat*FIRE_CO2_PER_ACRE,
                           CO2_No_Treat = Acres_Burned_No_Treat*FIRE_CO2_PER_ACRE,
                           CO2_Realized = Acres_Burned_Realized*FIRE_CO2_PER_ACRE)

#### Health Costs

treated_data_imp <- mutate(treated_data_imp, 
                           Deaths_With_Treat = Acres_Burned_With_Treat*Deaths_PER_ACRE_pred,
                           Deaths_No_Treat = Acres_Burned_No_Treat*Deaths_PER_ACRE_pred,
                           Deaths_Realized = Acres_Burned_Realized*Deaths_PER_ACRE_pred,
                           DeathCost_With_Treat = Acres_Burned_With_Treat*DeathCost_PER_ACRE_pred,
                           DeathCost_No_Treat = Acres_Burned_No_Treat*DeathCost_PER_ACRE_pred,
                           DeathCost_Realized = Acres_Burned_Realized*DeathCost_PER_ACRE_pred,
                           EarningsLoss_With_Treat = Acres_Burned_With_Treat*EarningsLoss_PER_ACRE_pred,
                           EarningsLoss_No_Treat = Acres_Burned_No_Treat*EarningsLoss_PER_ACRE_pred,
                           EarningsLoss_Realized = Acres_Burned_Realized*EarningsLoss_PER_ACRE_pred)

Death_savings <- (sum(treated_data_imp$DeathCost_No_Treat, na.rm = T) - sum(treated_data_imp$DeathCost_With_Treat, na.rm = T)) # Deaths*VSL
Death_savings_pct <- (sum(treated_data_imp$DeathCost_No_Treat, na.rm = T) - sum(treated_data_imp$DeathCost_With_Treat, na.rm = T))/sum(treated_data_imp$DeathCost_Realized, na.rm = T)

Death_savings


Earnings_savings <- (sum(treated_data_imp$EarningsLoss_No_Treat, na.rm = T) - sum(treated_data_imp$EarningsLoss_With_Treat, na.rm = T)) # Earnings Estimate
Earnings_savings_pct <- (sum(treated_data_imp$EarningsLoss_No_Treat, na.rm = T) - sum(treated_data_imp$EarningsLoss_With_Treat, na.rm = T))/sum(treated_data_imp$EarningsLoss_Realized, na.rm = T)

Earnings_savings

Health_Savings <- Death_savings + Earnings_savings

Avoided_Damages_imp <- house_val_savings + carbon_savings + Health_Savings

treated_data_imp <- mutate(treated_data_imp, lambda = case_when(TREAT_SIZE < q1 ~ lambda_1,
                                                                TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ lambda_2,
                                                                TREAT_SIZE > q2 ~ lambda_3,))

treated_data_imp <- mutate(treated_data_imp, Size_Bin = case_when(TREAT_SIZE < q1 ~ "Small",
                                                                  TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ "Medium",
                                                                  TREAT_SIZE > q2 ~ "Large",))

treat_ids <- unique(treated_data_imp$TREAT_ID)
treat_ids <- treat_ids[!is.na(treat_ids)]

treatment_costs <- treated_data_imp %>%
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
Avoided_Damages_NoImp <- house_val_savings + carbon_savings + Health_Savings
Benefit_Cost_Ratio_Conditional_NoImp = Avoided_Damages_NoImp/Treat_Cost

Treat_Size_Benefits_imp <- treated_data_imp %>%
  group_by(Size_Bin) %>%
  summarise(Acres_Burned_Savings = sum(Acres_Burned_No_Treat - Acres_Burned_With_Treat),
            Struc_Savings = sum(Structure_Burned_No_Treat - Structure_Burned_With_Treat),
            House_Savings = sum(Structure_Value_Lost_No_Treat - Structure_Value_Lost_With_Treat, na.rm = T),
            Carbon_Savings = (sum(CO2_No_Treat - CO2_With_Treat))*185,
            Death_Savings = sum(DeathCost_No_Treat - DeathCost_With_Treat, na.rm = T),
            Earning_Savings = sum(EarningsLoss_No_Treat - EarningsLoss_With_Treat, na.rm = T),
            N_Treats = n_distinct(TREAT_ID),
            lambda = dplyr::first(lambda))

Treat_Size_Benefits_imp <- mutate(Treat_Size_Benefits_imp, Total_Benefit = House_Savings + Carbon_Savings + Death_Savings)
Treat_Size_Benefits_imp <- mutate(Treat_Size_Benefits_imp, Benefit_Per_Treat = Total_Benefit/N_Treats)

#### Benefit per Treatment

sum(Treat_Size_Benefits_imp$Total_Benefit)/sum(Treat_Size_Benefits_imp$N_Treats)

#### Benefit Cost Ratios for different sized treatments

mu_1_imp <- filter(Treat_Size_Benefits_imp, Size_Bin == "Small")$Benefit_Per_Treat
mu_2_imp <- filter(Treat_Size_Benefits_imp, Size_Bin == "Medium")$Benefit_Per_Treat
mu_3_imp <- filter(Treat_Size_Benefits_imp, Size_Bin == "Large")$Benefit_Per_Treat

#### Average Benefit Per Acres

Treat_Sizes <- treated_data_imp %>%
  group_by(TREAT_ID) %>%
  summarise(Size_Bin = dplyr::first(Size_Bin), 
            TREAT_SIZE = dplyr::first(TREAT_SIZE),
  )

Size_Bin_Acres_Treated <- Treat_Sizes %>%
  group_by(Size_Bin) %>%
  summarise(Total_Acres_Treated = sum(TREAT_SIZE))

Treat_Size_Benefits_imp <- merge(Treat_Size_Benefits_imp, Size_Bin_Acres_Treated, by = "Size_Bin", all.x = T)

Treat_Size_Benefits_imp <- mutate(Treat_Size_Benefits_imp, Benefit_Per_Acre = Total_Benefit/Total_Acres_Treated)

#### Calculate the Ex-ante B/C Ratio for treatments completed 2007-2023

facts_act_grouped_06_23 <- mutate(facts_act_grouped_06_23, lambda = case_when(TREAT_SIZE < q1 ~ lambda_1,
                                                                              TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ lambda_2,
                                                                              TREAT_SIZE > q2 ~ lambda_3),
                                  mu = case_when(TREAT_SIZE < q1 ~ mu_1_imp,
                                                 TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ mu_2_imp,
                                                 TREAT_SIZE > q2 ~ mu_3_imp))

facts_act_grouped_06_23 <- mutate(facts_act_grouped_06_23, Size_Bin = case_when(TREAT_SIZE < q1 ~ "Small",
                                                                                TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ "Medium",
                                                                                TREAT_SIZE > q2 ~ "Large"))


facts_act_grouped_06_23 <- mutate(facts_act_grouped_06_23, Benefits = mu*lambda, 
                                  Benefit_Cost = Benefits/TOT_COST)

# facts_act_grouped_06_23 <- filter(facts_act_grouped_06_23, ACCOMPLISHED_ACRES > 0)

## Identify observations with high cost/acre

mean_cost <- mean(facts_west$COST_PER_UOM, na.rm = TRUE)
sd_cost <- sd(facts_west$COST_PER_UOM, na.rm = TRUE)
threshold <- mean_cost + 10*sd_cost


facts_act_grouped_06_23_filt <- filter(facts_act_grouped_06_23, COST_PER_UOM <= threshold)
facts_act_grouped_06_23_filt_1 <- filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q_05)

BCR_Overall_imp <- sum(facts_act_grouped_06_23_filt_1$Benefits)/sum(facts_act_grouped_06_23_filt_1$TOT_COST)

#### Ex ante B-C Ratio by Treatment Size Bins: "Small" (75-600), "Medium" (600-2400), "Large" (>2400)

N_Treats_imp <- sum(Treat_Size_Benefits_imp$N_Treats)

Robustness[5, ] <- c(
  round(BCR_Overall_imp, digits = 2),
  N_Treats_imp,
  "2017-2020"
)

write_csv(Robustness, here("data", "temp", "BCR_Robustness.csv"))






# BC_Large_noimp <- sum(filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q2)$Benefits)/sum(filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q2)$TOT_COST)
# BC_Medium_noimp <- sum(filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q1 & TREAT_SIZE <= q2)$Benefits)/sum(filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q1 & TREAT_SIZE <= q2)$TOT_COST)
# BC_Small_noimp <- sum(filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q_05 & TREAT_SIZE <= q1)$Benefits)/sum(filter(facts_act_grouped_06_23_filt, TREAT_SIZE > q_05 & TREAT_SIZE <= q1)$TOT_COST)
# 
# c(BC_Large_noimp, BC_Medium_noimp, BC_Small_noimp)
# 
# Robustness[3, ] <- c(
#   round(BCR_Overall_noimp, digits = 2),
#   round(BC_Small_noimp, digits = 2),
#   round(BC_Medium_noimp, digits = 2), 
#   round(BC_Large_noimp, digits = 2), 
#   paste0("$", prettyNum(round(  Avoided_Damages_Noimp, digits = 2), scientific = F, big.mark = ",")),
#   N_Treats_noimp
# )
# 
# write_csv(Robustness, here("data", "temp", "BCR_Robustness.csv"))





# #### Distribution for mu for each treatment size
# 
# Treatment_Benefits <- treated_data_filt %>%
#   group_by(TREAT_ID) %>%
#   summarise(House_Savings = sum(Structure_Value_Lost_No_Treat - Structure_Value_Lost_With_Treat, na.rm = T),
#             Carbon_Savings = (sum(CO2_No_Treat - CO2_With_Treat))*185,
#             Death_Savings = sum(DeathCost_No_Treat - DeathCost_With_Treat, na.rm = T),
#             Earning_Savings = sum(EarningsLoss_No_Treat - EarningsLoss_With_Treat, na.rm = T),
#             Size_Bin = dplyr::first(Size_Bin),
#             TREAT_SIZE = dplyr::first(TREAT_SIZE))
# 
# Treatment_Benefits <- mutate(Treatment_Benefits, Benefits = House_Savings + Carbon_Savings + Death_Savings + Earning_Savings)
# 
# Treatment_Benefits <- mutate(Treatment_Benefits, Benefit_Per_Acre = Benefits/TREAT_SIZE)
# 
# Treatment_Benefits$Size_Bin <- factor(Treatment_Benefits$Size_Bin, levels = c("Small", "Medium", "Large"))
# 
# Treatment_Benefits <- mutate(Treatment_Benefits, Benefits_rep = ifelse(Benefits < 0, 0, Benefits),
#                              Log_Benefit = sign(Benefit_Per_Acre) * log(abs(Benefit_Per_Acre) + 1))
# 
# Treatment_Benefits_plot <- ggplot(Treatment_Benefits, aes(x = Benefit_Per_Acre, fill = Size_Bin)) +
#   geom_density(alpha = 0.6) +  # Semi-transparent density curves
#   scale_fill_manual(
#     values = nature_palette, 
#     labels = c("Small (75-600)", "Medium (600-2400)", "Large (> 2400)")
#   ) +
#   labs(
#     x = "Log Treatment Benefit Per Acre Treated ($)", 
#     y = "Density",
#     fill = NULL  # Remove default legend title
#   ) +
#   theme_minimal(base_size = 16) +
#   theme(
#     panel.grid = element_blank(),  # Remove grid lines
#     axis.title = element_text(size = 18, face = "bold"),  
#     axis.text = element_text(size = 15),  
#     legend.position = "bottom",
#     legend.text = element_text(size = 15),
#     legend.key = element_rect(fill = NA),
#     legend.spacing.x = unit(0.3, 'cm'),
#     plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
#     panel.border = element_rect(color = "black", fill = NA, size = 1)
#   )
# 
# Treatment_Benefits_plot
# 
# 
# 
# 
# 
# #### Distribution of B-C Ratio based on our Treatment Size Bins: "Small" (75-600), "Medium" (600-2400), "Large" (>2400)
# 
# facts_act_grouped_17_23_filt <- mutate(facts_act_grouped_17_23_filt, BC = Benefits/TOT_COST)
# 
# facts_act_grouped_17_23_filt_1 <- filter(facts_act_grouped_17_23_filt, TREAT_SIZE > q_05)
# 
# facts_act_grouped_17_23_filt_1$Size_Bin <- factor(facts_act_grouped_17_23_filt_1$Size_Bin, levels = c("Small", "Medium", "Large"))
# 
# Median_BCR <- median(facts_act_grouped_17_23_filt_1$Benefit_Cost)
# Total_BCR <- sum(facts_act_grouped_17_23_filt_1$Benefits)/sum(facts_act_grouped_17_23_filt_1$TOT_COST)
# 
# 
# # Calculate median and total BCR
# Median_BCR <- median(facts_act_grouped_17_23_filt_1$Benefit_Cost, na.rm = TRUE)
# Total_BCR <- sum(facts_act_grouped_17_23_filt_1$Benefits, na.rm = TRUE) / sum(facts_act_grouped_17_23_filt_1$TOT_COST, na.rm = TRUE)
# 
# # Subset data to avoid dropping Size_Bin
# filtered_data <- facts_act_grouped_17_23_filt_1 %>% filter(BC < 25)
# 
# vline_data <- data.frame(
#   xintercept = c(1, Median_BCR, Total_BCR),
#   label = c("Cost-Effective", "Median BCR", "Total BCR"),
#   color = c("black", "#9C2A2F", "blue")  # Black, Gold, Deep Blue
# )
# 
# # Calculate the maximum density value for the BC variable
# max_density <- max(density(filter(facts_act_grouped_17_23_filt_1, BC < 25)$BC)$y)
# 
# # Create the plot
# Ex_Ante_BCR_plot <- ggplot(filter(facts_act_grouped_17_23_filt_1, BC < 25), aes(x = BC, fill = Size_Bin)) +
#   geom_density(alpha = 0.6) +  # Semi-transparent density curves
#   
#   # Vertical lines with explicit data, avoiding inheritance issues
#   geom_vline(data = vline_data, aes(xintercept = xintercept, color = color), 
#              linetype = "dashed", size = 1, inherit.aes = FALSE) +  
#   
#   # Labels for vertical lines, adjusting position closer to the x-axis
#   geom_text(data = vline_data, aes(x = xintercept, y = max_density * 0.10, label = label, color = color), 
#             angle = 90, vjust = -0.5, hjust = 0, size = 5, fontface = "bold", inherit.aes = FALSE) +
#   
#   scale_fill_manual(
#     values = nature_palette, 
#     labels = c("Small (75-600)", "Medium (600-2400)", "Large (> 2400)")
#   ) +
#   scale_color_identity() +  # Ensures colors from vline_data are used correctly
#   labs(
#     x = "Ex-Ante B-C Ratio", 
#     y = "Density",
#     fill = NULL  # Remove default legend title
#   ) +
#   theme_minimal(base_size = 16) +
#   theme(
#     panel.grid = element_blank(),  # Remove grid lines
#     axis.title = element_text(size = 18, face = "bold"),  
#     axis.text = element_text(size = 15),  
#     legend.position = "bottom",
#     legend.text = element_text(size = 15),
#     legend.key = element_rect(fill = NA),
#     legend.spacing.x = unit(0.3, 'cm'),
#     plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
#     panel.border = element_rect(color = "black", fill = NA, size = 1)
#   ) +
#   coord_cartesian(clip = "off")  # Prevents cropping of axis lines
# 
# Ex_Ante_BCR_plot
# 
# 
# #### Create small table of Ex-Ante B-C's
# 
# BC_df <- data.frame(
#   BC_1 = BC_1,
#   BC_2 = BC_2,
#   BC_3 = BC_3,
#   BC_Ratio_Ex_Ante_17_23 = BC_Ratio_Ex_Ante_17_23,
#   Median_BC = Median_BC
# )
# 
# colnames(BC_df) <- c("Small (75-600)","Medium (600-2400)", "Large (>2400)", "Total", "Median")
# 
# rownames(BC_df) <- "Ex-Ante Benefit Cost Ratios"
# 
# BC_df[1,] <- round(BC_df[1,], digits = 2)
# 
# table <- tableGrob(BC_df)

# ##### Combine all three plots into one
# 
# # Theme for table (same as before)
# clean_theme <- ttheme_minimal(
#   core = list(
#     fg_params = list(fontsize = 12),
#     bg_params = list(fill = NA, col = "grey80")
#   ),
#   colhead = list(
#     fg_params = list(fontsize = 12, fontface = "bold"),
#     bg_params = list(fill = NA, col = "grey80")
#   ),
#   padding = unit(c(6, 6), "mm")
# )
# 
# # Table object
# table <- tableGrob(BC_df, rows = NULL, theme = clean_theme)
# 
# # Panel tag and title in a horizontal layout
# panel_label <- textGrob("d)", x = 0, hjust = 0, gp = gpar(fontsize = 12, fontface = "bold"))
# table_title <- textGrob("Summary of Ex-Ante Benefit Cost Ratios", x = 0, hjust = 0, gp = gpar(fontsize = 13, fontface = "bold"))
# 
# # Combine d) and title in a row
# title_row <- arrangeGrob(
#   panel_label, table_title, ncol = 2, widths = unit.c(unit(0.05, "npc"), unit(0.95, "npc"))
# )
# 
# # Combine title row and table in a vertical layout
# table_with_title <- arrangeGrob(
#   title_row,
#   table,
#   ncol = 1,
#   heights = unit.c(unit(0.2, "npc"), unit(1, "npc"))
# )

# Final layout, balanced heights
# final_layout <- (
#   (lambda_plot | Treatment_Benefits_plot) / 
#     Ex_Ante_BCR_plot) +  
#   plot_annotation(
#     theme = theme(
#       plot.tag = element_text(face = "bold", size = 12),
#       plot.title = element_text(size = 16, face = "bold")
#     )
#   )
# 
# final_layout
# 
# ggsave(here("output", "figures", "Figure5.pdf"), final_layout, width = 14, height = 10, units = "in")







# pdf(here("output", "figures", "lambda_size_plot.pdf"))
# 
# print(lambda_plot)
# 
# dev.off()
# 
# pdf(here("output", "figures", "Treatment_Benefits_plot.pdf"))
# 
# print(Treatment_Benefits_plot)
# 
# dev.off()
# 
# pdf(here("output", "figures", "Ex_Ante_BCR_plot.pdf"))
# 
# print(Ex_Ante_BCR_plot)
# 
# dev.off()




#### Descriptive Stats of fires

HealthStats <- Grids_df_filt_1 %>%
  group_by(YEAR) %>%
  summarise(
    DeathCosts     = sum(DeathCost) / 1e12,
    EarningsLosses = sum(EarningsLoss) / 1e12,
    PM25_YEARLY = sum(FIRE_PM25)
  ) %>%
  mutate(
    HealthLosses = DeathCosts + EarningsLosses
  )

sum(filter(HealthStats, YEAR <= 2020)$HealthLosses)
sum(filter(HealthStats, YEAR > 2020)$HealthLosses)

sum(filter(HealthStats, YEAR <= 2020)$PM25_YEARLY)
sum(filter(HealthStats, YEAR > 2020)$PM25_YEARLY)



