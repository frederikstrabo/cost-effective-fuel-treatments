##########################################    02_conditional_effects.R ##########################################

################# Purpose: Run baseline spatial DiD regressions, create event study plots, and heterogeneity analysis.

################# Outputs: Figure 3 saved as "Figure3.pdf" saved in "output/figures".
#################          Table S1 saved as "TableS1.tex" saved in "output/tables".
#################          Figure S7 saved as "FigureS7.pdf" saved in "output/figures".
#################          Figure 5 saved as "Figure5.pdf" saved in "output/figures".
#################          Figure S5 saved as "FigureS5.pdf" saved in "output/figures".

################# Note: event study coefficient estimates might slightly change depending on which version of didimputation and fixest installed.

################# Estimated run time: ~2 minutes


rm(list=ls())

# if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr,sf, tmap, magrittr, rnaturalearth, rnaturalearthdata, ggplot2, maps, lwgeom, raster, stars, haven, stargazer, quantmod, lubridate, tidyr, ggpubr, 
               exactextractr, tictoc, terra, gtools, here, fixest, modelsummary, readr, rdrobust, prism, parallel,tmaptools, 
               OpenStreetMap, maptiles, gifski, geosphere, did, cowplot, didimputation, patchwork, stringr)

# Set Path
here::i_am("code/06_analysis/02_conditional_effects.R")

#### Load plot panel  

Grids_df <- read_csv(here("data", "intermediate", "SpatialDiD_Grids_L24_K05.csv"))
Grids_df <- mutate(Grids_df,  dist_treated = ifelse(Treated_Dir == 0, 1000, L_TAU)) # For Sun & Abraham Estimator
# Define treatment and control observations: post-treatment observations are considered treated if they are inside of a treatment (Treated) or directly after one (Treated_LAG) 
Grids_df_treated <- filter(Grids_df, time_to_treat > 0) %>% filter(Treated == 1 | Treated_LAG == 1) 
Grids_df_yet_treated <- filter(Grids_df, time_to_treat <= 0)

Grids_df <- rbind(Grids_df_treated, Grids_df_yet_treated)

###### Run Event Study DiD Regressions

# subset to observations that are "yet-to-be extinguished" in order to estimate the conditional hazard model
Grids_df_filt <- filter(Grids_df, BURN_LAG == 1) 

# limit to the 2.5km window - also remove treated observations that occur in distance bin 1 (i.e. Burn_Right_Away == 1).

Grids_df_filt_1 <- filter(Grids_df_filt, time_to_treat %in% seq(-5, 4,1) & Burn_Right_Away == 0)
Grids_df_filt_1_BURN <- filter(Grids_df_filt_1, BURN == 1) # Conditional on burning sample used to estimate burn severity effects

# Define controls

controls <- c("Slope", "Elev", "TRI", "Distance_FS_Road", "USFS_NF", "Wilderness", "MFRI", 
              "Distance_US_Highway", "Distance_WUI", "ERC", "WindSpeed", "FM1000", "PREV_BURN_10Y",
              "LAT", "LAT_LINE_INT", "LOG_DIST_LAT", "WIND_DIFF", "DELTA_MTT", "NA_DELTA_MTT",
              "LOG_FIRE_INTENSITY", "ROAD", "WUI")

# Controls in heterogeneity analysis - including fireline controls

controls_full <- c("Slope", "Elev", "TRI", "Distance_FS_Road", "USFS_NF", "Wilderness", "MFRI", 
                   "Distance_US_Highway", "Distance_WUI", "ERC", "WindSpeed", "FM1000", "PREV_BURN_10Y",
                   "LAT", "LAT_LINE_INT", "LOG_DIST_LAT", "WIND_DIFF", "DELTA_MTT", "NA_DELTA_MTT",
                   "LOG_FIRE_INTENSITY", "ROAD", "WUI",  "SUP_LINE", "SUP_LINE_INT", "LOG_SUP_LINE_DIST")

# Define fixed effects

FEs <- c("distance_bin", "direction_fire_FE", "FuelType_2001", "Aspect_Class")

# Define formula used 

first_stage_formula <- paste(
  paste(controls, collapse = " + "),  # Controls
  "|",  
  paste(FEs, collapse = " + ")  # Fixed Effects
)

#### Example of how to run event study in using did_imputation package - i.e. Borusyak et al. (2024) imputation method ####

es <- did_imputation(
  data = Grids_df_filt_1, yname = "BURN", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE,
  pretrends = -4:-1
)

# print(es, n = 100)

es_filt <- es %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

# Calculate the number of observations
N_obs_bin <- Grids_df_filt_1 %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Merge with the coefficient data
es_filt <- es_filt %>%
  left_join(N_obs_bin, by = "term")

interval_labels <- c(
  "-2.5 to -2", "-2 to -1.5", "-1.5 to -1", "-1 to -0.5",
  "-0.5 to 0", "0 to 0.5", "0.5 to 1", "1 to 1.5",
  "1.5 to 2", "2 to 2.5"
)

ggplot(es_filt, aes(x = term, y = estimate)) +
  # Main points
  geom_point(size = 3, color = "#0072B2") +
  # Error bars
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.15, linewidth = 0.8, color = "#0072B2") +
  # Dashed reference line at 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  # Dashed vertical line at -0.5
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "black", linewidth = 0.8) +
  # Observation count above each point
  geom_text(aes(y = conf.high + 0.02, label = N_obs),
            size = 4, color = "black", fontface = "bold") +

  # Axis styling
  scale_x_continuous(breaks = seq(min(es_filt$term), max(es_filt$term), 1),
                     labels = interval_labels,
                     name = "Distance from Treatment (km)") +
  scale_y_continuous(name = "Coefficient Estimate") +
  # Minimal but structured theme with axes
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),  # Remove minor gridlines for a cleaner look
    panel.border = element_rect(color = "black", fill = NA, linewidth = .4),  # Add axes
    axis.line = element_blank(),  # Avoid duplicate axis lines
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.ticks = element_line(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12)
  )


# #### Aside what is R^2 with and without MTT outputs ####
# 
# Grids_df_filt_1_test <- Grids_df_filt_1
# 
# Grids_df_filt_1_test$LOG_FIRE_INTENSITY[Grids_df_filt_1_test$LOG_FIRE_INTENSITY == -Inf] <- NA
# 
# Grids_df_filt_1_test <- mutate(Grids_df_filt_1_test, LOG_FIRE_INTENSITY_SCALE = scale(LOG_FIRE_INTENSITY), 
#                                DELTA_MTT_SCALE = scale(DELTA_MTT))
# 
# Grids_df_filt_1_test_BURN <- filter(Grids_df_filt_1_test, BURN == 1)
# 
# controls_scaled <- c("Slope", "Elev", "TRI", "Distance_FS_Road", "USFS_NF", "Wilderness", "MFRI", 
#                      "Distance_US_Highway", "Distance_WUI", "ERC", "WindSpeed", "FM1000", "PREV_BURN_10Y",
#                      "LAT", "LAT_LINE_INT", "LOG_DIST_LAT", "WIND_DIFF", "DELTA_MTT_SCALE", "NA_DELTA_MTT",
#                      "LOG_FIRE_INTENSITY_SCALE", "ROAD", "WUI",  "SUP_LINE", "SUP_LINE_INT", "LOG_SUP_LINE_DIST")
# 
# first_stage_formula_scaled <- paste(
#   paste(controls_scaled, collapse = " + "),  # Controls
#   "|",  
#   paste(FEs, collapse = " + ")  # Fixed Effects
# )
# 
# es_test <- did_imputation(
#   data = Grids_df_filt_1_test, yname = "BURN", gname = "L_TAU",
#   tname = "distance_bin", idname = "direction_fire_FE",
#   cluster_var = "FIRE_ID",
#   first_stage = as.formula(paste("~", first_stage_formula_scaled)),  # Use as.formula
#   wname = "Grid_Acres",
#   # event-study
#   horizon = TRUE, 
#   pretrends = -4:-1
# )
# 
# print(es_test, n = 100)
# 
# 
# es_test$p_value <- 2 * (1 - pnorm(abs(es_test$estimate / es_test$std.error)))
# 
# es_test %>% 
#   filter(term %in% c("DELTA_MTT_SCALE", "LOG_FIRE_INTENSITY_SCALE", "NA_DELTA_MTT")) %>%
#   mutate(
#     p_value_fmt = format(p_value, scientific = FALSE, digits = 4)
#   ) %>%
#   dplyr::select(term, estimate, std.error, p_value_fmt)
# 
# es_test_BS <- did_imputation(
#   data = Grids_df_filt_1_test, yname = "BURN_SEV_MED_HIGH_PCT", gname = "L_TAU",
#   tname = "distance_bin", idname = "direction_fire_FE",
#   cluster_var = "FIRE_ID",
#   first_stage = as.formula(paste("~", first_stage_formula_scaled)),  # Use as.formula
#   wname = "Grid_Acres",
#   # event-study
#   horizon = TRUE, 
#   pretrends = -4:-1
# )
# 
# print(es_test_BS, n = 100)
# 
# es_test_BS$p_value <- 2 * (1 - pnorm(abs(es_test_BS$estimate / es_test_BS$std.error)))
# 
# es_test_BS %>% 
#   filter(term %in% c("DELTA_MTT_SCALE", "LOG_FIRE_INTENSITY_SCALE", "NA_DELTA_MTT")) %>%
#   mutate(
#     p_value_fmt = format(p_value, scientific = FALSE, digits = 4)
#   ) %>%
#   dplyr::select(term, estimate, std.error, p_value_fmt)
# 


################ Run Event Study Plots and create Figure 3 ################

#### Function "Event_Study_Plot" takes as input: the sample used, dependent variable, its name, and the Xkm window to produce an event study plot figure 

Event_Study_Plot <- function(sample, dep_var, dep_var_name, window, OBS) {
  
  min_window <- min(window)
  
  y_axis_label <- "Treatment Effect"
  
  es <- did_imputation(
    data = sample, yname = dep_var, gname = "L_TAU",
    tname = "distance_bin", idname = "direction_fire_FE",
    cluster_var = "FIRE_ID",
    first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
    wname = "Grid_Acres",
    # event-study
    horizon = TRUE,
    pretrends = min_window:-1
  )
  
  # print(es, n = 100)
  
  es_filt <- es %>%
    filter(term %in% window) %>%
    mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
    bind_rows(tibble(lhs = "BURN", term = paste(min_window - 1), estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
    mutate(term = as.numeric(term))  # Convert to numeric for correct ordering
  
  # Calculate the number of observations
  N_obs_bin <- sample %>%
    filter(Treated_Dir == 1) %>%
    group_by(time_to_treat) %>%
    summarise(N_obs = n(), .groups = "drop") %>%
    mutate(term = time_to_treat)
  
  # Merge with the coefficient data
  es_filt <- es_filt %>%
    left_join(N_obs_bin, by = "term")
  
  n_obs_lab = es_filt$conf.high + 0.02
  
  if (dep_var == "BURN" | dep_var == "BURN_SEV_MED_HIGH_PCT" | dep_var == "BURN_SEV_HIGH_PCT"){
    
    es_filt$estimate <-  es_filt$estimate*100
    es_filt$std.error <-  es_filt$std.error*100
    es_filt$conf.low <-  es_filt$conf.low*100
    es_filt$conf.high <-  es_filt$conf.high*100
    
    y_axis_label <- "Treatment Effect (% points)"
    
    n_obs_lab = es_filt$conf.high + 1
  }
  
  # interval_labels <- c(
  #   "-2.5 to -2", "-2 to -1.5", "-1.5 to -1", "-1 to -0.5",
  #   "-0.5 to 0", "0 to 0.5", "0.5 to 1", "1 to 1.5",
  #   "1.5 to 2", "2 to 2.5"
  # )
  
  interval_labels <- c(
    "-2.5 to\n -2",
    "-2 to\n -1.5",
    "-1.5 to\n -1",
    "-1 to\n -0.5",
    "-0.5 to\n 0",
    "0 to\n 0.5",
    "0.5 to\n 1",
    "1 to\n 1.5",
    "1.5 to\n 2",
    "2 to\n 2.5"
  )
  
  if (length(window) == 15){
    
    interval_labels <- c(
      "-4 to\n -3.5",
      "-3.5 to\n -3",
      "-3 to\n -2.5",
      "-2.5 to\n -2",
      "-2 to\n -1.5",
      "-1.5 to\n -1",
      "-1 to\n -0.5",
      "-0.5 to\n 0",
      "0 to\n 0.5",
      "0.5 to\n 1",
      "1 to\n 1.5",
      "1.5 to\n 2",
      "2 to\n 2.5",
      "2.5 to\n 3",
      "3 to\n 3.5",
      "3.5 to\n 4"
    )
    
  }
  
  if (length(window) == 29){
    
    interval_labels <- c(
      "-7.5 to\n -7",
      "-7 to\n -6.5",
      "-6.5 to\n -6",
      "-6 to\n -5.5",
      "-5.5 to\n -5",
      "-5 to\n -4.5",
      "-4.5 to\n -4",
      "-4 to\n -3.5",
      "-3.5 to\n -3",
      "-3 to\n -2.5",
      "-2.5 to\n -2",
      "-2 to\n -1.5",
      "-1.5 to\n -1",
      "-1 to\n -0.5",
      "-0.5 to\n 0",
      "0 to\n 0.5",
      "0.5 to\n 1",
      "1 to\n 1.5",
      "1.5 to\n 2",
      "2 to\n 2.5",
      "2.5 to\n 3",
      "3 to\n 3.5",
      "3.5 to\n 4",
      "4 to\n 4.5",
      "4.5 to\n 5",
      "5 to\n 5.5",
      "5.5 to\n 6",
      "6 to\n 6.5",
      "6.5 to\n 7",
      "7 to\n 7.5"
    )
    
  }
  
  if (OBS == F){
    
    # Create the plot
    plot <- ggplot(es_filt, aes(x = term, y = estimate)) +
      # Main points
      geom_point(size = 3, color = "#0072B2") +
      # Error bars
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                    width = 0.15, linewidth = 0.8, color = "#0072B2") +
      # Dashed reference line at 0
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
      # Dashed vertical line at -0.5
      geom_vline(xintercept = -0.5, linetype = "dashed", color = "black", linewidth = 0.8) +
      # Observation count above each point
      # geom_text(aes(y = conf.high + 0.02, label = N_obs),
      #           size = 4, color = "black", fontface = "bold") +
      
      # Axis styling
      scale_x_continuous(breaks = seq(min(es_filt$term), max(es_filt$term), 1),
                         labels = interval_labels,
                         name = "Distance from treatment interaction (km)") +
      scale_y_continuous(name = y_axis_label) +
      
      # Minimal but structured theme with axes
      theme_minimal(base_family = "Helvetica") +
      theme(
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
        axis.title = element_text(size = 8, face = "bold"),
        axis.text  = element_text(size = 7),
        plot.title   = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                    margin = margin(b = 2)),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(3, 3, 3, 3)
      ) +
      ggtitle(dep_var_name)
    
  }
  
  if (OBS == T){
    
    # Create the plot
    plot <- ggplot(es_filt, aes(x = term, y = estimate)) +
      # Main points
      geom_point(size = 3, color = "#0072B2") +
      # Error bars
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                    width = 0.15, linewidth = 0.8, color = "#0072B2") +
      # Dashed reference line at 0
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
      # Dashed vertical line at -0.5
      geom_vline(xintercept = -0.5, linetype = "dashed", color = "black", linewidth = 0.8) +
      # Observation count above each point
      geom_text(aes(y = n_obs_lab, label = N_obs),
                size = 4, color = "black", fontface = "bold") +
      
      # Axis styling
      scale_x_continuous(breaks = seq(min(es_filt$term), max(es_filt$term), 1),
                         labels = interval_labels,
                         name = "Distance from treatment interaction (km)") +
      scale_y_continuous(name = y_axis_label) +
      
      # Minimal but structured theme with axes
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
      )+
      ggtitle(dep_var_name)
    
  }
  
  
  
  return(plot)
}

# #### Alternative Event study plot that uses grey-scale coloring
# 
# Event_Study_Plot <- function(sample, dep_var, dep_var_name, window, OBS) {
#   
#   min_window <- min(window)
#   y_axis_label <- "Treatment Effect"
#   
#   es <- did_imputation(
#     data = sample, yname = dep_var, gname = "L_TAU",
#     tname = "distance_bin", idname = "direction_fire_FE",
#     cluster_var = "FIRE_ID",
#     first_stage = as.formula(paste("~", first_stage_formula)),
#     wname = "Grid_Acres",
#     horizon = TRUE,
#     pretrends = min_window:-1
#   )
#   
#   es_filt <- es %>%
#     dplyr::filter(term %in% window) %>%
#     dplyr::mutate(term = as.character(term)) %>%
#     dplyr::bind_rows(tibble::tibble(lhs = "BURN", term = paste(min_window - 1),
#                                     estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
#     dplyr::mutate(term = as.numeric(term))
#   
#   N_obs_bin <- sample %>%
#     dplyr::filter(Treated_Dir == 1) %>%
#     dplyr::group_by(time_to_treat) %>%
#     dplyr::summarise(N_obs = dplyr::n(), .groups = "drop") %>%
#     dplyr::mutate(term = time_to_treat)
#   
#   es_filt <- es_filt %>%
#     dplyr::left_join(N_obs_bin, by = "term")
#   
#   n_obs_lab <- es_filt$conf.high + 0.02
#   
#   if (dep_var %in% c("BURN", "BURN_SEV_MED_HIGH_PCT", "BURN_SEV_HIGH_PCT")) {
#     es_filt$estimate  <- es_filt$estimate  * 100
#     es_filt$std.error <- es_filt$std.error * 100
#     es_filt$conf.low  <- es_filt$conf.low  * 100
#     es_filt$conf.high <- es_filt$conf.high * 100
#     
#     y_axis_label <- "Treatment Effect (% points)"
#     n_obs_lab <- es_filt$conf.high + 1
#   }
#   
#   # --- your interval_labels code unchanged ---
#   interval_labels <- c(
#     "-2.5 to\n -2",
#     "-2 to\n -1.5",
#     "-1.5 to\n -1",
#     "-1 to\n -0.5",
#     "-0.5 to\n 0",
#     "0 to\n 0.5",
#     "0.5 to\n 1",
#     "1 to\n 1.5",
#     "1.5 to\n 2",
#     "2 to\n 2.5"
#   )
#   if (length(window) == 15) { interval_labels <- c(
#     "-4 to\n -3.5","-3.5 to\n -3","-3 to\n -2.5","-2.5 to\n -2","-2 to\n -1.5",
#     "-1.5 to\n -1","-1 to\n -0.5","-0.5 to\n 0","0 to\n 0.5","0.5 to\n 1",
#     "1 to\n 1.5","1.5 to\n 2","2 to\n 2.5","2.5 to\n 3","3 to\n 3.5","3.5 to\n 4"
#   )}
#   if (length(window) == 29) { interval_labels <- c(
#     "-7.5 to\n -7","-7 to\n -6.5","-6.5 to\n -6","-6 to\n -5.5","-5.5 to\n -5",
#     "-5 to\n -4.5","-4.5 to\n -4","-4 to\n -3.5","-3.5 to\n -3","-3 to\n -2.5",
#     "-2.5 to\n -2","-2 to\n -1.5","-1.5 to\n -1","-1 to\n -0.5","-0.5 to\n 0",
#     "0 to\n 0.5","0.5 to\n 1","1 to\n 1.5","1.5 to\n 2","2 to\n 2.5","2.5 to\n 3",
#     "3 to\n 3.5","3.5 to\n 4","4 to\n 4.5","4.5 to\n 5","5 to\n 5.5","5.5 to\n 6",
#     "6 to\n 6.5","6.5 to\n 7","7 to\n 7.5"
#   )}
#   
#   # --- grayscale-safe styling ---
#   ink       <- "grey20"  # main series (points + CI)
#   ref_ink   <- "black"   # reference lines
#   ci_lwd    <- 0.8
#   pt_size   <- 3
#   
#   base_plot <- ggplot2::ggplot(es_filt, ggplot2::aes(x = term, y = estimate)) +
#     ggplot2::geom_errorbar(
#       ggplot2::aes(ymin = conf.low, ymax = conf.high),
#       width = 0.15, linewidth = ci_lwd, color = ink
#     ) +
#     ggplot2::geom_point(size = pt_size, color = ink) +
#     ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = ref_ink, linewidth = 0.6) +
#     ggplot2::geom_vline(xintercept = -0.5, linetype = "dashed", color = ref_ink, linewidth = 0.6) +
#     ggplot2::scale_x_continuous(
#       breaks = seq(min(es_filt$term), max(es_filt$term), 1),
#       labels = interval_labels,
#       name = "Distance from treatment interaction (km)"
#     ) +
#     ggplot2::scale_y_continuous(name = y_axis_label) +
#     ggplot2::theme_minimal(base_family = "Helvetica") +
#     ggplot2::theme(
#       panel.grid = ggplot2::element_blank(),
#       panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.4),
#       axis.title = ggplot2::element_text(size = 8, face = "bold"),
#       axis.text  = ggplot2::element_text(size = 7),
#       plot.title = ggplot2::element_text(size = 8.5, face = "bold", hjust = 0.5,
#                                          margin = ggplot2::margin(b = 2)),
#       axis.ticks = ggplot2::element_line(linewidth = 0.3),
#       plot.margin = ggplot2::margin(3, 3, 3, 3)
#     ) +
#     ggplot2::ggtitle(dep_var_name)
#   
#   if (OBS) {
#     base_plot <- base_plot +
#       ggplot2::geom_text(ggplot2::aes(y = n_obs_lab, label = N_obs),
#                          size = 4, color = "black", fontface = "bold") +
#       ggplot2::theme(
#         axis.text.x = ggplot2::element_text(
#           angle = 45, hjust = 1, vjust = 1, size = 6, lineheight = 0.9
#         )
#       )
#   }
#   
#   base_plot
# }



#### Create event study plots in our baseline sample using different outcomes

plot1 <- Event_Study_Plot(Grids_df_filt_1, "BURN", "Conditional Burn Probability", seq(-4,4,1), F)
plot2 <- Event_Study_Plot(Grids_df_filt_1_BURN, "BURN_SEV0", "Conditional Burn Severity", seq(-4,4,1), F)
plot3 <- Event_Study_Plot(Grids_df_filt_1_BURN, "BURN_SEV_MED_HIGH_PCT", "Conditional Moderate–High Severity Burn", seq(-4,4,1), F)
plot4 <- Event_Study_Plot(Grids_df_filt_1_BURN, "BURN_SEV_HIGH_PCT", "Conditional High Severity Burn", seq(-4,4,1), F)

saveRDS(plot1, here("output", "event_plotBURN.rds"))
saveRDS(plot2, here("output", "event_plotBURNSEV.rds"))
saveRDS(plot3, here("output", "event_plotBURNSEV_MED_HIGH.rds"))
saveRDS(plot4, here("output", "event_plotBURNSEV_HIGH.rds"))


#### Format Figure 3

# Combine the two plots and ensure only one legend
final_layout <- (plot1 | plot3) +
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(
      family = "Helvetica",
      face   = "bold",
      size   = 10
    )
  ) &
  theme(
    panel.spacing = unit(0.8, "lines")   # tighter spacing
  )

ggsave(here("output", "figures", "Figure3.pdf"),
       final_layout,
       width  = 7.24,          # Science 3-column width
       height = 4.5,    # ~4.22 in
       units  = "in",
       dpi    = 300)

################ Create Baseline Table of Results - Table S1 ################

## run_fixest_burn is a function that will run the model using fixest - necessary for initializing results in table used by the etable package

run_fixest_burn <- function(data) {
  feols(BURN ~ sunab(dist_treated, distance_bin) + 
          Slope + Elev + TRI + factor(Aspect_Class) + factor(FuelType_2001) + Distance_FS_Road + MFRI + 
          Distance_US_Highway + Distance_WUI + ERC + WindSpeed + FM1000 + PREV_BURN_10Y + LAT + LAT_LINE_INT + 
          LOG_DIST_LAT + WIND_DIFF + DELTA_MTT + NA_DELTA_MTT + LOG_FIRE_INTENSITY + ROAD + WUI + Wilderness | 
          direction_fire_FE + distance_fire_FE,
        cluster = ~FIRE_ID,  
        weights = ~Grid_Acres,
        data = data)
}

## Function that will take estimates from the borusyak model and convert them into the fixest model for table outputting

update_fixest_coefs <- function(fixest_model, es_filt) {
  new_coefs <- es_filt %>%
    filter(term %in% 0:4) %>%
    arrange(term)
  
  coef_names <- paste0("distance_bin::", new_coefs$term)
  
  new_coef_table <- coeftable(fixest_model)
  
  new_coef_table[coef_names, "Estimate"] <- new_coefs$estimate
  new_coef_table[coef_names, "Std. Error"] <- new_coefs$std.error
  new_coef_table[coef_names, "Pr(>|t|)"] <- 2 * pt(abs(new_coefs$estimate / new_coefs$std.error), 
                                                   df = df.residual(fixest_model), 
                                                   lower.tail = FALSE)
  
  fixest_model$coeftable <- new_coef_table
  return(fixest_model)
}

event_twfe_burn <- run_fixest_burn(Grids_df_filt_1)
event_twfe_burn <- update_fixest_coefs(event_twfe_burn, es_filt)


## Separate function for burn severity

run_fixest_burnsev <- function(data) {
  feols(BURN_SEV_MED_HIGH_PCT ~ sunab(dist_treated, distance_bin) + 
          Slope + Elev + TRI + factor(Aspect_Class) + factor(FuelType_2001) + Distance_FS_Road + MFRI + 
          Distance_US_Highway + Distance_WUI + ERC + WindSpeed + FM1000 + PREV_BURN_10Y + LAT + LAT_LINE_INT + 
          LOG_DIST_LAT + WIND_DIFF + DELTA_MTT + NA_DELTA_MTT + LOG_FIRE_INTENSITY + ROAD + WUI + Wilderness | 
          direction_fire_FE + distance_fire_FE,
        cluster = ~FIRE_ID,  
        weights = ~Grid_Acres,
        data = data)
}

run_did_imputation <- function(data, var) {
  es <- did_imputation(
    data = data, yname = var, gname = "L_TAU",
    tname = "distance_bin", idname = "direction_fire_FE",
    cluster_var = "FIRE_ID",
    first_stage = as.formula(paste("~", first_stage_formula)),
    wname = "Grid_Acres",
    horizon = TRUE, 
    pretrends = -4:-1
  )
  
  es_filt <- es %>%
    filter(term %in% seq(-4, 4, 1)) %>%
    mutate(term = as.character(term)) %>%
    bind_rows(tibble(lhs = "BURN", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
    mutate(term = as.numeric(term))
  
  return(es_filt)
}

event_twfe_burnsev <- run_fixest_burnsev(Grids_df_filt_1_BURN)
imputation_results_burnsev <- run_did_imputation(Grids_df_filt_1_BURN, "BURN_SEV_MED_HIGH_PCT")
event_twfe_burnsev <- update_fixest_coefs(event_twfe_burnsev, imputation_results_burnsev)

imputation_results_High_burnsev <- run_did_imputation(Grids_df_filt_1_BURN, "BURN_SEV_HIGH_PCT")
event_twfe_High_burnsev <- update_fixest_coefs(event_twfe_burnsev, imputation_results_High_burnsev)

## Define variables for etable

var_dict <- c(
  BURN = "Probability of Fire Spread",
  BURN_SEV0 = "Conditional Burn Severity",
  BURN_SEV_MED_HIGH_PCT = "Conditional Moderate-High Burn Severity",
  direction_fire_FE = "Direction-Fire",
  distance_fire_FE = "Distance-Fire",
  `distance_bin::0` = "$Treat_0$",
  `distance_bin::1` = "$Treat_1$",
  `distance_bin::2` = "$Treat_2$",
  `distance_bin::3` = "$Treat_3$",
  `distance_bin::4` = "$Treat_4$",
  `distance_bin::5` = "$Treat_5$"
)

setFixest_dict(var_dict)

spatDiDmods <- list()
spatDiDmods[["Burn Probability - (1)"]] <- event_twfe_burn
spatDiDmods[["Burn Probability - (2)"]] <- event_twfe_burnsev


## Create & save table of baseline results for Table S1

baseline_results <- etable(spatDiDmods,
                           keep = c("%distance_bin::0", "%distance_bin::1", "%distance_bin::2", "%distance_bin::3", "%distance_bin::4"),  # Use original names with "%"
                           digits = "r3",
                           digits.stats = 2,
                           fitstat = c("n", "r2"),
                           style.tex = style.tex("aer",
                                                 fixef.suffix = " FEs",  # Remove the " FEs" suffix
                                                 fixef.where = "var", # Keep fixed effects where they are placed
                                                 yesNo = c("Yes", "No")
                           ),
                           signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.1), # Custom stars
                           tex = TRUE)
baseline_results

baseline_results %>% write_lines(here("output", "tables", "TableS1.tex"))

## p-values reported in the paper

reg_summary <- summary(event_twfe_burn)
p_values <- reg_summary$coeftable[, "Pr(>|t|)"]
round(p_values, digits = 5)

reg_summary_BS <- summary(event_twfe_burnsev)
p_values_BS <- reg_summary_BS$coeftable[, "Pr(>|t|)"]
round(p_values_BS, digits = 5)

## Percentage effects

mean_BS <- mean(filter(Grids_df_filt_1_BURN, Treated == 0)$BURN_SEV_MED_HIGH_PCT, na.rm = T)
reg_summary_BS$coeftable[c("distance_bin::0", "distance_bin::1", "distance_bin::2", "distance_bin::3", "distance_bin::4"), "Estimate"]/mean_BS


reg_summary_burn <- summary(event_twfe_burn)
reg_summary_burn$coeftable[c("distance_bin::0", "distance_bin::1", "distance_bin::2", "distance_bin::3", "distance_bin::4"), "Estimate"]

# Percent High Burn Severity Impacts

reg_summary_High_BS <- summary(event_twfe_High_burnsev)
mean_High_BS <- mean(filter(Grids_df_filt_1_BURN, Treated == 0)$BURN_SEV_HIGH_PCT, na.rm = T)
reg_summary_High_BS$coeftable[c("distance_bin::0", "distance_bin::1", "distance_bin::2", "distance_bin::3", "distance_bin::4"), "Estimate"]/mean_High_BS



################    Estimate results with different distance windows - Figure S7   ################

#### Different Event Windows
Grids_df_filt_1_4km <- filter(Grids_df_filt, time_to_treat %in% seq(-8, 7,1) & Burn_Right_Away == 0)
Grids_df_filt_1_7km <- filter(Grids_df_filt, time_to_treat %in% seq(-15, 14,1) & Burn_Right_Away == 0)

plot1 <- Event_Study_Plot(Grids_df_filt_1, "BURN", "Burn Probability", seq(-4,4,1), F)
plot2 <- Event_Study_Plot(Grids_df_filt_1_4km, "BURN", "Burn Probability", seq(-7,7,1), F)
plot3 <- Event_Study_Plot(Grids_df_filt_1_7km, "BURN", "Burn Probability", seq(-14,14,1), F)

#### Distribution of observations based on event times

Grids_df_Treated <- filter(Grids_df_filt, Burn_Right_Away == 0 & Treated_Dir == 1) %>% filter(time_to_treat %in% seq(-15, 14,1))

dist_time_to_treat_plot <- ggplot(Grids_df_Treated, aes(x = time_to_treat)) +
  # Histogram with refined colors
  geom_histogram(binwidth = 1, color = "gray30", fill = "#377EB8", alpha = 0.8) +
  
  # Data labels with better alignment
  stat_bin(
    binwidth = 1,
    geom = "text",
    aes(label = after_stat(count),
        y = after_stat(count) + max(after_stat(count)) * 0.03),
    size = 2.5,      # smaller text
    fontface = "bold",
    color = "black",
    hjust = 0.5
  ) +
  scale_x_continuous(
    breaks = seq(-15, 14, 1),
    labels = seq(-15, 14, 1),
    name   = "Distance from treatment interaction (0.5 km bins)"
  ) +
  
  # Titles and labels
  labs(
    x = "Distance from treatment interaction (0.5 km bins)",
    y = "Number of observations"
  ) +
  
  # Professional theme with refinements
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

#### Combine Plots

plot1 <- plot1 +
  scale_x_continuous(
    breaks = seq(-5,4,1),
    labels = seq(-5,4,1),
    name = "Distance from treatment interaction (0.5 km bins)"
  )

plot2 <- plot2 +
  scale_x_continuous(
    breaks = seq(-8,7,1),
    labels = seq(-8,7,1),
    name = "Distance from treatment interaction (0.5 km bins)"
  )

plot3 <- plot3 +
  scale_x_continuous(
    breaks = seq(-15,14,1),
    labels = seq(-15,14,1),
    name = "Distance from treatment interaction (0.5 km bins)"
  )


## Layout of Figure S7

final_layout <- dist_time_to_treat_plot /
  (plot1 | plot2) /
  plot3 +
  plot_layout(
    heights = c(1.3, 1, 1),   # A taller than rows 2 & 3
    guides  = "collect"
  ) +
  plot_annotation(tag_levels = "A") &   # adds A, B, C, D tags
  theme(
    plot.tag = element_text(
      family = "Helvetica",
      face   = "bold",
      size   = 10
    ),
    panel.spacing = unit(0.6, "lines"),
    legend.position = "bottom",
    legend.title = element_text(family = "Helvetica", face = "bold", size = 7),
    legend.text  = element_text(family = "Helvetica", size = 7)
  )


ggsave(here("output", "figures", "FigureS7.pdf"), final_layout,        
       width  = 7.24,          # Science 3-column width
       height = 6,    # ~4.5 in
       units  = "in",
       dpi    = 300)



####################################     Heterogeneity Analysis  - Create Figure 5    ####################################  

###### Heterogeneity By Suppression Effort

## Define an observation as close to suppression effort if it's within 0.2 km of a LAT drop or suppression line 

Grids_df_filt_1 <- mutate(Grids_df_filt_1, EFFORT = ifelse(DIST_LAT <= 0.20 | SUP_LINE_DIST <= 0.20, 1, 0))

## Filter to fires with fire line info

Grids_df_filt_1_lines <- filter(Grids_df_filt_1, SUP_LINES_PRESENT == 1)

#### Define samples where all the treated observations are close to suppression effort and unrestricted controls - we call this "Grids_df_filt_1_SE"

Grids_df_filt_1_effort <- filter(Grids_df_filt_1_lines, Treated == 1 | Treated_LAG == 1) %>% filter(EFFORT == 1)
Grids_df_filt_1_effort_BURN <- filter(Grids_df_filt_1_effort, BURN == 1)
Grids_df_filt_1_controls <- filter(Grids_df_filt_1_lines, Treated == 0 & Treated_LAG == 0) 

Grids_df_filt_1_SE <- rbind(Grids_df_filt_1_effort, Grids_df_filt_1_controls)
Grids_df_filt_1_SE_BURN <- filter(Grids_df_filt_1_SE, BURN == 1)


#### Define other sample where all the treated observations are not close to suppression effort and unrestricted controls - we call this "Grids_df_filt_1_NOSE"

Grids_df_filt_1_noeffort <- filter(Grids_df_filt_1_lines, Treated == 1 | Treated_LAG == 1) %>% filter(EFFORT == 0) 
Grids_df_filt_1_noeffort_BURN <- filter(Grids_df_filt_1_noeffort, BURN == 1)
Grids_df_filt_1_controls <- filter(Grids_df_filt_1_lines, Treated == 0 & Treated_LAG == 0) 

Grids_df_filt_1_NOSE <- rbind(Grids_df_filt_1_noeffort, Grids_df_filt_1_controls)
Grids_df_filt_1_NOSE_BURN <- filter(Grids_df_filt_1_NOSE, BURN == 1)


#### Estimate DiD model under the different samples separately

first_stage_formula_suppression <- paste(
  paste(controls_full, collapse = " + "),  # Controls
  "|",  
  paste(FEs, collapse = " + ")  # Fixed Effects
)

es_sup_effort <- did_imputation(
  data = Grids_df_filt_1_SE, yname = "BURN", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula_suppression)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_sup_effort_BS <- did_imputation(
  data = Grids_df_filt_1_SE_BURN, yname = "BURN_SEV0", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula_suppression)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_sup_noeffort <- did_imputation(
  data = Grids_df_filt_1_NOSE, yname = "BURN", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula_suppression)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_sup_noeffort_BS <- did_imputation(
  data = Grids_df_filt_1_NOSE_BURN, yname = "BURN_SEV_MED_HIGH_PCT", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula_suppression)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_sup_effort_filt <- es_sup_effort %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

es_sup_effort_filt_BS <- es_sup_effort_BS %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN_SEV_MED_HIGH_PCT", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

es_sup_noeffort_filt <- es_sup_noeffort %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

es_sup_noeffort_filt_BS <- es_sup_noeffort_BS %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN_SEV_MED_HIGH_PCT", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering


# Calculate the number of observations in each event time
N_obs_bin_effort <- Grids_df_filt_1_effort %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

N_obs_bin_effort_BS <- Grids_df_filt_1_effort_BURN %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Calculate the number of observations in each event time
N_obs_bin_noeffort <- Grids_df_filt_1_noeffort %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Calculate the number of observations
N_obs_bin_noeffort_BS <- Grids_df_filt_1_noeffort_BURN %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)


# Combine datasets and add group labels
es_combined <- bind_rows(
  es_sup_effort_filt %>% mutate(Group = "With Suppression"),
  es_sup_noeffort_filt %>% mutate(Group = "No Suppression")
)

es_combined_BS <- bind_rows(
  es_sup_effort_filt_BS %>% mutate(Group = "With Suppression"),
  es_sup_noeffort_filt_BS %>% mutate(Group = "No Suppression")
)

# Merge with observation count
N_obs_combined <- bind_rows(
  N_obs_bin_effort %>% mutate(Group = "With Suppression"),
  N_obs_bin_noeffort %>% mutate(Group = "No Suppression")
)

# Merge with observation count
N_obs_combined_BS <- bind_rows(
  N_obs_bin_effort_BS %>% mutate(Group = "With Suppression"),
  N_obs_bin_noeffort_BS %>% mutate(Group = "No Suppression")
)

es_combined <- es_combined %>%
  left_join(N_obs_combined, by = c("term", "Group"))

es_combined_BS <- es_combined_BS %>%
  left_join(N_obs_combined_BS, by = c("term", "Group"))

# Define colors
colors <- c("With Suppression" = "#0072B2", "No Suppression" = "#D55E00")


es_combined <- es_combined %>%
  mutate(Group = factor(Group, levels = c("With Suppression", "No Suppression")),
         term_offset = ifelse(Group == "With Suppression", term - 0.2, term + 0.2)
  )

es_combined_BS <- es_combined_BS %>%
  mutate(Group = factor(Group, levels = c("With Suppression", "No Suppression")),
         term_offset = ifelse(Group == "With Suppression", term - 0.2, term + 0.2)
  )

es_combined <- filter(es_combined, term %in% seq(0,4,1))
es_combined_BS <- filter(es_combined_BS, term %in% seq(0,4,1))

interval_labels <- c(
  "0 to 0.5", "0.5 to 1", "1 to 1.5",
  "1.5 to 2", "2 to 2.5"
)

# Create the plot
sup_plot <- ggplot(es_combined, aes(x = term_offset, y = estimate*100, color = Group, group = Group)) +
  # Main points
  geom_point(size = 3) +
  # Error bars
  geom_errorbar(aes(ymin = conf.low*100, ymax = conf.high*100), 
                width = 0.15, linewidth = 0.8) +
  # Error bars
  geom_errorbar(aes(ymin = conf.low*100, ymax = conf.high*100), width = 0.2, linewidth = 0.8) +
  # Dashed reference line at 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Observation count above each point
  # geom_text(aes(y = conf.high*100 + 0.02, label = N_obs), size = 3.5, color = "black") +
  
  # Axis styling
  scale_x_continuous(breaks = seq(min(es_combined$term), max(es_combined$term), 1),
                     name = "Distance from treatment interaction (km)", 
                     labels = interval_labels) +
  scale_y_continuous(name = "Treatment Effect (% points)") +
  scale_color_manual(name = "Suppression Effort", values = colors) +  # Change legend title
  # Add title
  ggtitle("Suppression Effort") +  # Minimal but structured theme with axes
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text  = element_text(size = 7),
    plot.title   = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                margin = margin(b = 2)),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(3, 3, 3, 3)
  )

sup_plot

# Create the plot
sup_plot_BS <- ggplot(es_combined_BS, aes(x = term_offset, y = estimate*100, color = Group, group = Group)) +
  # Main points
  geom_point(size = 3) +
  # Error bars
  geom_errorbar(aes(ymin = conf.low*100, ymax = conf.high*100), 
                width = 0.15, linewidth = 0.8) +
  # Error bars
  geom_errorbar(aes(ymin = conf.low*100, ymax = conf.high*100), width = 0.2, linewidth = 0.8) +
  # Dashed reference line at 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Observation count above each point
  # geom_text(aes(y = conf.high + 0.02, label = N_obs), size = 3.5, color = "black") +
  
  # Axis styling
  scale_x_continuous(breaks = seq(min(es_combined$term), max(es_combined$term), 1),
                     name = "Distance from treatment interaction (km)", 
                     labels = interval_labels) +
  scale_y_continuous(name = "Treatment Effect (% points)") +
  scale_color_manual(name = "Suppression Effort", values = colors) +  # Change legend title
  # Add title
  ggtitle("Suppression Effort") +  # Minimal but structured theme with axes
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text  = element_text(size = 7),
    plot.title   = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                margin = margin(b = 2)),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(3, 3, 3, 3)
  )

sup_plot_BS



###### Heterogeneity by treatment type


Grids_df_filt_1_Mech <- filter(Grids_df_filt_1, Treated == 1 | Treated_LAG == 1) %>% filter(TREAT_CAT == "Mechanical")
Grids_df_filt_1_controls <- filter(Grids_df_filt_1, Treated == 0 & Treated_LAG == 0)

Grids_df_filt_1_Mech <- rbind(Grids_df_filt_1_Mech, Grids_df_filt_1_controls)
Grids_df_filt_1_Mech_BURN <- filter(Grids_df_filt_1_Mech, BURN == 1)

Grids_df_filt_1_RX <- filter(Grids_df_filt_1, Treated == 1 | Treated_LAG == 1) %>% filter(TREAT_CAT == "Fire")
Grids_df_filt_1_controls <- filter(Grids_df_filt_1, Treated == 0 & Treated_LAG == 0)

Grids_df_filt_1_RX <- rbind(Grids_df_filt_1_RX, Grids_df_filt_1_controls)
Grids_df_filt_1_RX_BURN <- filter(Grids_df_filt_1_RX, BURN == 1)

Grids_df_filt_1_RX_Mech <- filter(Grids_df_filt_1, Treated == 1 | Treated_LAG == 1) %>% filter(TREAT_CAT == "Rx & Mech")
Grids_df_filt_1_controls <- filter(Grids_df_filt_1, Treated == 0 & Treated_LAG == 0)

Grids_df_filt_1_RX_Mech <- rbind(Grids_df_filt_1_RX_Mech, Grids_df_filt_1_controls)
Grids_df_filt_1_RX_Mech_BURN <- filter(Grids_df_filt_1_RX_Mech, BURN == 1)


es_mech <- did_imputation(
  data = Grids_df_filt_1_Mech, yname = "BURN", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_mech_BS <- did_imputation(
  data = Grids_df_filt_1_Mech_BURN, yname = "BURN_SEV_MED_HIGH_PCT", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_rx <- did_imputation(
  data = Grids_df_filt_1_RX, yname = "BURN", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_rx_BS <- did_imputation(
  data = Grids_df_filt_1_RX_BURN, yname = "BURN_SEV_MED_HIGH_PCT", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_rx_mech <- did_imputation(
  data = Grids_df_filt_1_RX_Mech, yname = "BURN", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_rx_mech_BS <- did_imputation(
  data = Grids_df_filt_1_RX_Mech_BURN, yname = "BURN_SEV_MED_HIGH_PCT", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_mech_filt <- es_mech %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

es_mech_BS_filt <- es_mech_BS %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN_SEV_MED_HIGH_PCT", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering


es_rx_filt <- es_rx %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

es_rx_BS_filt <- es_rx_BS %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN_SEV_MED_HIGH_PCT", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

es_rx_mech_filt <- es_rx_mech %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

es_rx_mech_BS_filt <- es_rx_mech_BS %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN_SEV_MED_HIGH_PCT", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

# Calculate the number of observations
N_obs_bin_mech <- Grids_df_filt_1_Mech %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Calculate the number of observations
N_obs_bin_mech_BS <- Grids_df_filt_1_Mech_BURN %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Calculate the number of observations
N_obs_bin_rx <- Grids_df_filt_1_RX %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Calculate the number of observations
N_obs_bin_rx_BS <- Grids_df_filt_1_RX_BURN %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Calculate the number of observations
N_obs_bin_rx_mech <- Grids_df_filt_1_RX_Mech %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Calculate the number of observations
N_obs_bin_rx_mech_BS <- Grids_df_filt_1_RX_Mech_BURN %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Combine datasets and add group labels
es_combined <- bind_rows(
  es_mech_filt %>% mutate(Group = "Mechanical Only"),
  es_rx_filt %>% mutate(Group = "RX Only"),
  es_rx_mech_filt %>% mutate(Group = "RX & Mech.")
)

# Combine datasets and add group labels
es_combined_BS <- bind_rows(
  es_mech_BS_filt %>% mutate(Group = "Mechanical Only"),
  es_rx_BS_filt %>% mutate(Group = "RX Only"),
  es_rx_mech_BS_filt %>% mutate(Group = "RX & Mech.")
)

# Merge with observation count
N_obs_combined <- bind_rows(
  N_obs_bin_mech %>% mutate(Group = "Mechanical Only"),
  N_obs_bin_rx %>% mutate(Group = "RX Only"),
  N_obs_bin_rx_mech %>% mutate(Group = "RX & Mech.")
)

# Merge with observation count
N_obs_combined_BS <- bind_rows(
  N_obs_bin_mech_BS %>% mutate(Group = "Mechanical Only"),
  N_obs_bin_rx_BS %>% mutate(Group = "RX Only"),
  N_obs_bin_rx_mech_BS %>% mutate(Group = "RX & Mech.")
)

es_combined <- es_combined %>%
  left_join(N_obs_combined, by = c("term", "Group"))

es_combined_BS <- es_combined_BS %>%
  left_join(N_obs_combined_BS, by = c("term", "Group"))

# Define colors
colors <- c("Mechanical Only" = "#0072B2", "RX Only" = "#D55E00", "RX & Mech." = "#33a02c")

# Offset terms to avoid overlap
es_combined <- es_combined %>%
  mutate(Group = factor(Group, levels = c("Mechanical Only", "RX Only", "RX & Mech.")),
         term_offset = case_when(
           Group == "Mechanical Only" ~ term - 0.2,
           Group == "RX Only" ~ term,
           Group == "RX & Mech." ~ term + 0.2
         ))

# Offset terms to avoid overlap
es_combined_BS <- es_combined_BS %>%
  mutate(Group = factor(Group, levels = c("Mechanical Only", "RX Only", "RX & Mech.")),
         term_offset = case_when(
           Group == "Mechanical Only" ~ term - 0.2,
           Group == "RX Only" ~ term,
           Group == "RX & Mech." ~ term + 0.2
         ))

es_combined <- filter(es_combined, term %in% seq(0,4,1))
es_combined_BS <- filter(es_combined_BS, term %in% seq(0,4,1))

# Create the plot
treat_type_plot <- ggplot(es_combined, aes(x = term_offset, y = estimate*100, color = Group, group = Group)) +
  # Main points
  geom_point(size = 3) +
  # Error bars
  geom_errorbar(aes(ymin = conf.low*100, ymax = conf.high*100), 
                width = 0.15, linewidth = 0.8) +
  # Error bars
  geom_errorbar(aes(ymin = conf.low*100, ymax = conf.high*100), width = 0.2, linewidth = 0.8) +
  # Dashed reference line at 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Observation count above each point
  # geom_text(aes(y = conf.high + 0.02, label = N_obs), size = 3.5, color = "black") +
  
  # Axis styling
  scale_x_continuous(breaks = seq(min(es_combined$term), max(es_combined$term), 1),
                     name = "Distance from treatment interaction (km)", 
                     labels = interval_labels) +
  scale_y_continuous(name = "Treatment Effect (% points)") +
  scale_color_manual(name = "Treatment Type", values = colors) +  # Change legend title
  # Add title
  ggtitle("Treatment Type") +  # Minimal but structured theme with axes
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text  = element_text(size = 7),
    plot.title   = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                margin = margin(b = 2)),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(3, 3, 3, 3)
  )

treat_type_plot


# Create the plot
treat_type_plot_BS <- ggplot(es_combined_BS, aes(x = term_offset, y = estimate*100, color = Group, group = Group)) +
  # Main points
  geom_point(size = 3) +
  # Error bars
  geom_errorbar(aes(ymin = conf.low*100, ymax = conf.high*100), 
                width = 0.15, linewidth = 0.8) +
  # Error bars
  geom_errorbar(aes(ymin = conf.low*100, ymax = conf.high*100), width = 0.2, linewidth = 0.8) +
  # Dashed reference line at 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Observation count above each point
  # geom_text(aes(y = conf.high + 0.02, label = N_obs), size = 3.5, color = "black") +
  
  # Axis styling
  scale_x_continuous(breaks = seq(min(es_combined$term), max(es_combined$term), 1),
                     name = "Distance from treatment interaction (km)", 
                     labels = interval_labels) +
  scale_y_continuous(name = "Treatment Effect (% points)") +
  scale_color_manual(name = "Treatment Type", values = colors) +  # Change legend title
  # Add title
  ggtitle("Treatment Type") +  # Minimal but structured theme with axes
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text  = element_text(size = 7),
    plot.title   = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                margin = margin(b = 2)),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(3, 3, 3, 3)
  )

treat_type_plot_BS


###### Heterogeneity by treatment size

quantile(Grids_df_filt_1$TREAT_SIZE, probs = c(.33, .66, 1), na.rm = T)


Grids_df_filt_1_small <- filter(Grids_df_filt_1, Treated == 1 | Treated_LAG == 1) %>% filter(TREAT_SIZE >= 75 & TREAT_SIZE <= 600)
Grids_df_filt_1_controls <- filter(Grids_df_filt_1, Treated == 0 & Treated_LAG == 0)

Grids_df_filt_1_small <- rbind(Grids_df_filt_1_small, Grids_df_filt_1_controls)
Grids_df_filt_1_small_BURN <- filter(Grids_df_filt_1_small, BURN == 1)


Grids_df_filt_1_medium <- filter(Grids_df_filt_1, Treated == 1 | Treated_LAG == 1) %>% filter(TREAT_SIZE > 600 & TREAT_SIZE <= 2400)
Grids_df_filt_1_controls <- filter(Grids_df_filt_1, Treated == 0 & Treated_LAG == 0)

Grids_df_filt_1_medium <- rbind(Grids_df_filt_1_medium, Grids_df_filt_1_controls)
Grids_df_filt_1_medium_BURN <- filter(Grids_df_filt_1_medium, BURN == 1)

Grids_df_filt_1_large <- filter(Grids_df_filt_1, Treated == 1 | Treated_LAG == 1) %>% filter(TREAT_SIZE > 2400)
Grids_df_filt_1_controls <- filter(Grids_df_filt_1, Treated == 0 & Treated_LAG == 0)

Grids_df_filt_1_large <- rbind(Grids_df_filt_1_large, Grids_df_filt_1_controls)
Grids_df_filt_1_large_BURN <- filter(Grids_df_filt_1_large, BURN == 1)


es_small <- did_imputation(
  data = Grids_df_filt_1_small, yname = "BURN", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_small_BS <- did_imputation(
  data = Grids_df_filt_1_small_BURN, yname = "BURN_SEV_MED_HIGH_PCT", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_medium <- did_imputation(
  data = Grids_df_filt_1_medium, yname = "BURN", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_medium_BS <- did_imputation(
  data = Grids_df_filt_1_medium_BURN, yname = "BURN_SEV_MED_HIGH_PCT", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_large <- did_imputation(
  data = Grids_df_filt_1_large, yname = "BURN", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_large_BS <- did_imputation(
  data = Grids_df_filt_1_large_BURN, yname = "BURN_SEV_MED_HIGH_PCT", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)


es_small_filt <- es_small %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

es_small_BS_filt <- es_small_BS %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN_SEV_MED_HIGH_PCT", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

es_medium_filt <- es_medium %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

es_medium_BS_filt <- es_medium_BS %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN_SEV_MED_HIGH_PCT", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

es_large_filt <- es_large %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

es_large_BS_filt <- es_large_BS %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN_SEV_MED_HIGH_PCT", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

# Calculate the number of observations
N_obs_bin_small <- Grids_df_filt_1_small %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Calculate the number of observations
N_obs_bin_small_BS <- Grids_df_filt_1_small_BURN %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Calculate the number of observations
N_obs_bin_medium <- Grids_df_filt_1_medium %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Calculate the number of observations
N_obs_bin_medium_BS <- Grids_df_filt_1_medium_BURN %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Calculate the number of observations
N_obs_bin_large <- Grids_df_filt_1_large %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Calculate the number of observations
N_obs_bin_large_BS <- Grids_df_filt_1_large_BURN %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Combine datasets and add group labels
es_combined <- bind_rows(
  es_small_filt %>% mutate(Group = "75-600 Acres"),
  es_medium_filt %>% mutate(Group = "600-2400 Acres"),
  es_large_filt %>% mutate(Group = "> 2400 Acres")
)

# Combine datasets and add group labels
es_combined_BS <- bind_rows(
  es_small_BS_filt %>% mutate(Group = "75-600 Acres"),
  es_medium_BS_filt %>% mutate(Group = "600-2400 Acres"),
  es_large_BS_filt %>% mutate(Group = "> 2400 Acres")
)

# Merge with observation count
N_obs_combined <- bind_rows(
  N_obs_bin_small %>% mutate(Group = "75-600 Acres"),
  N_obs_bin_medium %>% mutate(Group = "600-2400 Acres"),
  N_obs_bin_large %>% mutate(Group = "> 2400 Acres")
)

# Merge with observation count
N_obs_combined_BS <- bind_rows(
  N_obs_bin_small_BS %>% mutate(Group = "75-600 Acres"),
  N_obs_bin_medium_BS %>% mutate(Group = "600-2400 Acres"),
  N_obs_bin_large_BS %>% mutate(Group = "> 2400 Acres")
)

es_combined <- es_combined %>%
  left_join(N_obs_combined, by = c("term", "Group"))

es_combined_BS <- es_combined_BS %>%
  left_join(N_obs_combined_BS, by = c("term", "Group"))

# Define colors
colors <- c("75-600 Acres" = "#0072B2", "600-2400 Acres" = "#D55E00", "> 2400 Acres" = "#33a02c")

es_combined <- es_combined %>%
  mutate(Group = factor(Group, levels = c("75-600 Acres", "600-2400 Acres", "> 2400 Acres")),
         term_offset = case_when(
           Group == "75-600 Acres" ~ term - 0.2,
           Group == "600-2400 Acres" ~ term,
           Group == "> 2400 Acres" ~ term + 0.2
         ))

es_combined_BS <- es_combined_BS %>%
  mutate(Group = factor(Group, levels = c("75-600 Acres", "600-2400 Acres", "> 2400 Acres")),
         term_offset = case_when(
           Group == "75-600 Acres" ~ term - 0.2,
           Group == "600-2400 Acres" ~ term,
           Group == "> 2400 Acres" ~ term + 0.2
         ))

es_combined <- filter(es_combined, term %in% seq(0,4,1))
es_combined_BS <- filter(es_combined_BS, term %in% seq(0,4,1))

# Create the plot
treat_size_plot <- ggplot(es_combined, aes(x = term_offset, y = estimate*100, color = Group, group = Group)) +
  # Main points
  geom_point(size = 3) +
  # Error bars
  geom_errorbar(aes(ymin = conf.low*100, ymax = conf.high*100), 
                width = 0.15, linewidth = 0.8) +
  # Error bars
  geom_errorbar(aes(ymin = conf.low*100, ymax = conf.high*100), width = 0.2, linewidth = 0.8) +
  # Dashed reference line at 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Observation count above each point
  # geom_text(aes(y = conf.high + 0.02, label = N_obs), size = 3.5, color = "black") +
  
  # Axis styling
  scale_x_continuous(breaks = seq(min(es_combined$term), max(es_combined$term), 1),
                     name = "Distance from treatment interaction (km)", 
                     labels = interval_labels) +
  scale_y_continuous(name = "Treatment Effect (% points)") +
  scale_color_manual(name = "Treatment Size", values = colors) +  # Change legend title
  # Add title
  ggtitle("Treatment Size") +  # Minimal but structured theme with axes
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text  = element_text(size = 7),
    plot.title   = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                margin = margin(b = 2)),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(3, 3, 3, 3)
  )

treat_size_plot


# Create the plot
treat_size_plot_BS <- ggplot(es_combined_BS, aes(x = term_offset, y = estimate*100, color = Group, group = Group)) +
  # Main points
  geom_point(size = 3) +
  # Error bars
  geom_errorbar(aes(ymin = conf.low*100, ymax = conf.high*100), 
                width = 0.15, linewidth = 0.8) +
  # Error bars
  geom_errorbar(aes(ymin = conf.low*100, ymax = conf.high*100), width = 0.2, linewidth = 0.8) +
  # Dashed reference line at 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Observation count above each point
  # geom_text(aes(y = conf.high + 0.02, label = N_obs), size = 3.5, color = "black") +
  
  # Axis styling
  scale_x_continuous(breaks = seq(min(es_combined$term), max(es_combined$term), 1),
                     name = "Distance from treatment interaction (km)", 
                     labels = interval_labels) +
  scale_y_continuous(name = "Treatment Effect (% points)") +
  scale_color_manual(name = "Treatment Size", values = colors) +  # Change legend title
  # Add title
  ggtitle("Treatment Size") +  # Minimal but structured theme with axes
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text  = element_text(size = 7),
    plot.title   = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                margin = margin(b = 2)),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(3, 3, 3, 3)
  )

treat_size_plot_BS


###### Heterogeneity by Time Since Treatment

Grids_df_filt_1_young <- filter(Grids_df_filt_1, Treated == 1 | Treated_LAG == 1) %>% filter(YEARS_SINCE_TREAT <= 3)
Grids_df_filt_1_controls <- filter(Grids_df_filt_1, Treated == 0 & Treated_LAG == 0)

Grids_df_filt_1_young <- rbind(Grids_df_filt_1_young, Grids_df_filt_1_controls)
Grids_df_filt_1_young_BURN <- filter(Grids_df_filt_1_young, BURN == 1)

Grids_df_filt_1_old <- filter(Grids_df_filt_1, Treated == 1 | Treated_LAG == 1) %>% filter(YEARS_SINCE_TREAT > 3 & YEARS_SINCE_TREAT <= 7)
Grids_df_filt_1_old <- rbind(Grids_df_filt_1_old, Grids_df_filt_1_controls)
Grids_df_filt_1_old_BURN <- filter(Grids_df_filt_1_old, BURN == 1)

Grids_df_filt_1_oldest <- filter(Grids_df_filt_1, Treated == 1 | Treated_LAG == 1) %>% filter(YEARS_SINCE_TREAT > 7)
Grids_df_filt_1_oldest <- rbind(Grids_df_filt_1_oldest, Grids_df_filt_1_controls)
Grids_df_filt_1_oldest_BURN <- filter(Grids_df_filt_1_oldest, BURN == 1)

es_young <- did_imputation(
  data = Grids_df_filt_1_young, yname = "BURN", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_young_BS <- did_imputation(
  data = Grids_df_filt_1_young_BURN, yname = "BURN_SEV_MED_HIGH_PCT", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_old <- did_imputation(
  data = Grids_df_filt_1_old, yname = "BURN", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_old_BS <- did_imputation(
  data = Grids_df_filt_1_old_BURN, yname = "BURN_SEV_MED_HIGH_PCT", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_oldest <- did_imputation(
  data = Grids_df_filt_1_oldest, yname = "BURN", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_oldest_BS <- did_imputation(
  data = Grids_df_filt_1_oldest_BURN, yname = "BURN_SEV_MED_HIGH_PCT", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE, 
  pretrends = -4:-1
)

es_young_filt <- es_young %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

es_young_BS_filt <- es_young_BS %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN_SEV_MED_HIGH_PCT", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

es_old_filt <- es_old %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

es_old_BS_filt <- es_old_BS %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN_SEV_MED_HIGH_PCT", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

es_oldest_filt <- es_oldest %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

es_oldest_BS_filt <- es_oldest_BS %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN_SEV_MED_HIGH_PCT", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

# Calculate the number of observations
N_obs_bin_young <- Grids_df_filt_1_young %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Calculate the number of observations
N_obs_bin_young_BS <- Grids_df_filt_1_young_BURN %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Calculate the number of observations
N_obs_bin_old <- Grids_df_filt_1_old %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Calculate the number of observations
N_obs_bin_old_BS <- Grids_df_filt_1_old_BURN %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Calculate the number of observations
N_obs_bin_oldest <- Grids_df_filt_1_oldest %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Calculate the number of observations
N_obs_bin_oldest_BS <- Grids_df_filt_1_oldest_BURN %>%
  filter(Treated_Dir == 1) %>%
  group_by(time_to_treat) %>%
  summarise(N_obs = n(), .groups = "drop") %>%
  mutate(term = time_to_treat)

# Combine datasets and add group labels
es_combined <- bind_rows(
  es_young_filt %>% mutate(Group = "0-3 Years"),
  es_old_filt %>% mutate(Group = "4-7 Years"),
  es_oldest_filt %>% mutate(Group = "8-10 Years")
)

# Combine datasets and add group labels
es_combined_BS <- bind_rows(
  es_young_BS_filt %>% mutate(Group = "0-3 Years"),
  es_old_BS_filt %>% mutate(Group = "4-7 Years"),
  es_oldest_BS_filt %>% mutate(Group = "8-10 Years")
)

# Merge with observation count
N_obs_combined <- bind_rows(
  N_obs_bin_young %>% mutate(Group = "0-3 Years"),
  N_obs_bin_old %>% mutate(Group = "4-7 Years"),
  N_obs_bin_oldest %>% mutate(Group = "8-10 Years")
)

# Merge with observation count
N_obs_combined_BS <- bind_rows(
  N_obs_bin_young_BS %>% mutate(Group = "0-3 Years"),
  N_obs_bin_old_BS %>% mutate(Group = "4-7 Years"),
  N_obs_bin_oldest_BS %>% mutate(Group = "8-10 Years")
)

es_combined <- es_combined %>%
  left_join(N_obs_combined, by = c("term", "Group"))

es_combined_BS <- es_combined_BS %>%
  left_join(N_obs_combined, by = c("term", "Group"))

# Define colors
colors <- c("0-3 Years" = "#0072B2", "4-7 Years" = "#D55E00", "8-10 Years" = "#33a02c")

es_combined <- es_combined %>%
  mutate(Group = factor(Group, levels = c("0-3 Years", "4-7 Years", "8-10 Years")),
         term_offset = case_when(
           Group == "0-3 Years" ~ term - 0.2,
           Group == "4-7 Years" ~ term,
           Group == "8-10 Years" ~ term + 0.2
         ))

es_combined_BS <- es_combined_BS %>%
  mutate(Group = factor(Group, levels = c("0-3 Years", "4-7 Years", "8-10 Years")),
         term_offset = case_when(
           Group == "0-3 Years" ~ term - 0.2,
           Group == "4-7 Years" ~ term,
           Group == "8-10 Years" ~ term + 0.2
         ))

es_combined <- filter(es_combined, term %in% seq(0,4,1))
es_combined_BS <- filter(es_combined_BS, term %in% seq(0,4,1))

# Create the plot
treat_time_plot <- ggplot(es_combined, aes(x = term_offset, y = estimate*100, color = Group, group = Group)) +
  # Main points
  geom_point(size = 3) +
  # Error bars
  geom_errorbar(aes(ymin = conf.low*100, ymax = conf.high*100), 
                width = 0.15, linewidth = 0.8) +
  # Error bars
  geom_errorbar(aes(ymin = conf.low*100, ymax = conf.high*100), width = 0.2, linewidth = 0.8) +
  # Dashed reference line at 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Observation count above each point
  # geom_text(aes(y = conf.high + 0.02, label = N_obs), size = 3.5, color = "black") +
  
  # Axis styling
  scale_x_continuous(breaks = seq(min(es_combined$term), max(es_combined$term), 1),
                     name = "Distance from treatment interaction (km)", 
                     labels = interval_labels) +
  scale_y_continuous(name = "Treatment Effect (% points)") +
  scale_color_manual(name = "Treatment Time", values = colors) +  # Change legend title
  # Add title
  ggtitle("Treatment Time") +  # Minimal but structured theme with axes
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text  = element_text(size = 7),
    plot.title   = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                margin = margin(b = 2)),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(3, 3, 3, 3)
  )

treat_time_plot

# Create the plot
treat_time_plot_BS <- ggplot(es_combined_BS, aes(x = term_offset, y = estimate*100, color = Group, group = Group)) +
  # Main points
  geom_point(size = 3) +
  # Error bars
  geom_errorbar(aes(ymin = conf.low*100, ymax = conf.high*100), 
                width = 0.15, linewidth = 0.8) +
  # Error bars
  geom_errorbar(aes(ymin = conf.low*100, ymax = conf.high*100), width = 0.2, linewidth = 0.8) +
  # Dashed reference line at 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Observation count above each point
  # geom_text(aes(y = conf.high + 0.02, label = N_obs), size = 3.5, color = "black") +
  
  # Axis styling
  scale_x_continuous(breaks = seq(min(es_combined$term), max(es_combined$term), 1),
                     name = "Distance from treatment interaction (km)", 
                     labels = interval_labels) +
  scale_y_continuous(name = "Treatment Effect (% points)") +
  scale_color_manual(name = "Treatment Time", values = colors) +  # Change legend title
  # Add title
  ggtitle("Treatment Time") +  # Minimal but structured theme with axes
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text  = element_text(size = 7),
    plot.title   = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                margin = margin(b = 2)),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(3, 3, 3, 3)
  )

treat_time_plot_BS




########  Layout the four sources of heterogeneity into one plot to create Figure 5

final_layout <- (sup_plot | treat_size_plot)/(treat_type_plot + treat_time_plot) +   
  plot_annotation(tag_levels = 'A') & 
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
    legend.position = "bottom",
    plot.margin       = margin(t = 3.5, r = 17.5, b = 3.5, l = 3.5),  # more room on right
    legend.box.margin = margin(0, 10, 0, 0)                         # nudge legend left a bit
  )

ggsave(here("output", "figures", "Figure5.pdf"), final_layout,        
       width  = 7.24,          # Science 3-column width
       height = 7,    # ~4.5 in
       units  = "in",
       dpi    = 300)

########  Layout the four sources of burn severity heterogeneity into one plot to create Figure S5

final_layout <- (sup_plot_BS | treat_size_plot_BS)/(treat_type_plot_BS + treat_time_plot_BS) +   
  plot_annotation(tag_levels = 'A') & 
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
    legend.position = "bottom",
    plot.margin       = margin(t = 3.5, r = 17.5, b = 3.5, l = 3.5),  # more room on right
    legend.box.margin = margin(0, 10, 0, 0)                         # nudge legend left a bit
  )

ggsave(here("output", "figures", "FigureS5.pdf"), final_layout,     
       width  = 7.24,          # Science 3-column width
       height = 7,    # ~4.5 in
       units  = "in",
       dpi    = 300)