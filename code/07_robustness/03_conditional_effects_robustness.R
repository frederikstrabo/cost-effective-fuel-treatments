##########################################    03_conditional_effects_robustness.R ##########################################

################# Purpose: Run robustness checks on the conditional effects spatial DiD regressions reported in supplemental appendix.

################# Outputs: i) Figure S8 saved as "FigureS8.pdf" stored in "output/figures".
#################         ii) Figure S9 saved as "FigureS9.pdf" in "output/figures".
#################        iv) Table S3 saved as "TableS3.tex" in "output/tables".
#################        v) Table S2 saved as "TableS2.tex" in "output/tables".
#################        vi) Table S6 saved as "TableS6.tex" in "output/tables".
#################        vii) Table S5 saved as "TableS5.tex" in "output/tables".
#################        viii) Table S7 saved as "TableS7.tex" in "output/tables".
#################        ix) Table S4 saved as "TableS4.tex" in "output/tables".
#################        x) Table S8 saved as "TableS8.tex" in "output/tables".

################# Estimated run time: ~ 1.5 hours 

rm(list=ls())

# if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr,sf, tmap, magrittr, rnaturalearth, rnaturalearthdata, ggplot2, maps, lwgeom, raster, stars, haven, stargazer, quantmod, lubridate, tidyr, ggpubr, 
               exactextractr, tictoc, terra, gtools, here, fixest, modelsummary, readr, rdrobust, prism, parallel,tmaptools, 
               OpenStreetMap, maptiles, gifski, geosphere, cowplot, didimputation, Matching, cobalt, patchwork, rgenoud, did)

# Set Path
here::i_am("code/07_robustness/03_conditional_effects_robustness.R")

print(packageVersion("did"))

tic()

#########  Create Figure S8 Incomplete Treatments Robustness Check ###########

Grids_df <- read_csv(here("data", "intermediate", "SpatialDiD_Grids_24L_K05_Incomp.csv")) %>%
  arrange(direction_fire_FE, distance_bin)
Grids_df <- mutate(Grids_df, LOG_DIST_LAT = log(DIST_LAT + 1))

# Define treatment and control observations: post-treatment observations are considered treated if they are inside of a treatment (Treated) or directly after one (Treated_LAG)

Grids_df_treated <- filter(Grids_df, time_to_treat > 0) %>% filter(Treated == 1 | Treated_LAG == 1)
Grids_df_yet_treated <- filter(Grids_df, time_to_treat <= 0)

Grids_df <- rbind(Grids_df_treated, Grids_df_yet_treated)

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

FEs <- c("distance_bin", "direction_fire_FE", "FuelType_2001", "Aspect_Class")

first_stage_formula <- paste(
  paste(controls, collapse = " + "),  # Controls
  "|",
  paste(FEs, collapse = " + ")  # Fixed Effects
)


#### Example of how to run event study in using did_imputation package - i.e. Borusyak et al. (2024) imputation method

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
                     name = "Distance from Treatment (0.5 km)") +
  scale_y_continuous(name = "Coefficient Estimate") +

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
  )

#### Function "Event_Study_Plot" takes as input: the sample used, dependent variable, its name, and the Xkm window to produce an event study plot figure

Event_Study_Plot <- function(sample, dep_var, dep_var_name, window) {

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


  # if (dep_var == "BURN"){
  #
  #   es_filt$estimate <-  es_filt$estimate*100
  #   es_filt$std.error <-  es_filt$std.error*100
  #   es_filt$conf.low <-  es_filt$conf.low*100
  #   es_filt$conf.high <-  es_filt$conf.high*100
  #
  #   y_axis_label <- "Treatment Effect (%pt)"
  #
  #   n_obs_lab = es_filt$conf.high + 0.5
  #
  # }

  interval_labels <- c(
    "-2.5 to -2", "-2 to -1.5", "-1.5 to -1", "-1 to -0.5",
    "-0.5 to 0", "0 to 0.5", "0.5 to 1", "1 to 1.5",
    "1.5 to 2", "2 to 2.5"
  )

  # Create the plot
  plot <- ggplot(es_filt, aes(x = term, y = estimate*100)) +
    # Main points
    geom_point(size = 3, color = "#0072B2") +
    # Error bars
    geom_errorbar(aes(ymin = conf.low*100, ymax = conf.high*100),
                  width = 0.20, linewidth = 0.8, color = "#0072B2") +
    # Dashed reference line at 0
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    # Dashed vertical line at -0.5
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "black", linewidth = 0.8) +
    # Observation count above each point
    # geom_text(aes(y = n_obs_lab, label = N_obs),
    #           size = 4, color = "black", fontface = "bold") +

    # Axis styling
    scale_x_continuous(breaks = seq(min(es_filt$term), max(es_filt$term), 1),
                       labels = interval_labels,
                       name = "Distance from treatment interaction (km)") +
    scale_y_continuous(name = "Treatment Effect (% points)") +
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
      plot.margin = margin(3, 3, 3, 3)) +
    ggtitle(dep_var_name)

  return(plot)
}

plot1 <- Event_Study_Plot(Grids_df_filt_1, "BURN", "Conditional Burn Probability", seq(-4,4,1))
plot2 <- Event_Study_Plot(Grids_df_filt_1_BURN, "BURN_SEV_MED_HIGH_PCT", "Conditional Moderate-High Burn Severity", seq(-4,4,1))


## layout of figure S8

final_layout <- (plot1 | plot2) +
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

## Save figure S7

ggsave(here("output", "figures", "FigureS8.pdf"), final_layout,
       width  = 7.24,          # Science 3-column width
       height = 4.5,    # ~4.22 in
       units  = "in",
       dpi    = 300)


######### Create Figure S9 Matching Robustness Check ###########

###### Build Matched Subsample

Grids_df <- read_csv(here("data", "intermediate", "SpatialDiD_Grids_L24_K05.csv"))
Grids_df <- mutate(Grids_df,  dist_treated = ifelse(Treated_Dir == 0, 1000, L_TAU)) # For Sun & Abraham Estimator

# Define treatment and control observations: post-treatment observations are considered treated if they are inside of a treatment (Treated) or directly after one (Treated_LAG)

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

Grids_df_filt_1$distance_fire_FE <- as.numeric(as.factor(Grids_df_filt_1$distance_fire_FE))

## X is the matrix of controls that we'd like to create balance on

X  <- Grids_df_filt_1 %>% dplyr::select(distance_bin, USFS_NF, WindSpeed, ERC,
                                        DELTA_MTT, LOG_FIRE_INTENSITY, LOG_DIST_LAT, TRI)

Tr  <- Grids_df_filt_1$Treated_Dir == 1

Y  <- Grids_df_filt_1$BURN

## Matrix of controls to match on

BalanceMat <- cbind(Grids_df_filt_1$distance_bin, Grids_df_filt_1$USFS_NF, Grids_df_filt_1$WindSpeed, Grids_df_filt_1$ERC,
                    Grids_df_filt_1$DELTA_MTT, Grids_df_filt_1$LOG_FIRE_INTENSITY, Grids_df_filt_1$LOG_DIST_LAT, Grids_df_filt_1$TRI)

## Match observations using the GenMatch function from the "Matching" package

## Will match plots exactly based on their distance bin and whether or not it occurs in a National Forest or not. Plots will then inexactly matched to find the optimal
##  covariate balance across the most important determinants of fire spread: wind speed, ERC, Delta T, log fire intensity (from MTT),
##    log distance to nearest LAT drop, and topographic ruggedness index.

## Will use a population size of 100

gen1 <- GenMatch(Tr = Tr, X = X,
                 BalanceMatrix = BalanceMat,
                 estimand = "ATT",
                 M = 1,
                 replace=FALSE,
                 exact = c(T, T, F, F, F, F, F, F),
                 pop.size=100,
                 max.generations=20,
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


var_labels <- c(
  "USFS_NF" = "National Forest",
  "WindSpeed" = "Wind Speed",
  "ERC" = "ERC",
  "DELTA_MTT" = "Arrival Time",
  "LOG_FIRE_INTENSITY" = "Log Fire Intensity",
  "LOG_DIST_LAT" = "Log Distance LAT",
  "TRI" = "TRI"
)

## Create balance plot - Figure S9 A)

balance_plot <- love.plot(balance.table_tr,
                          threshold = .1,
                          line = TRUE,
                          stars = NULL,          # <- no stars
                          var.names = var_labels) +
  theme_minimal(base_size = 14, base_family = "serif") +  # Use a serif font
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +  # Orange & Blue (Nature-like)
  scale_linetype_manual(values = c("solid", "dashed")) +  # Make dashed threshold
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
  labs(
    title = "Covariate balance before and after matching",
    x = "Standardized mean difference (treated - untreated)",
    y = "Covariates",
    color = "Sample"
  )

## Estimate first using fixest to initialize regression output that will be used in etable - will replace coefficients with borusyak estimates

event_twfe_Matching = feols(BURN ~ sunab(dist_treated, distance_bin) + ## Our key interaction: time × treatment status
                              Slope + Elev + TRI + factor(Aspect_Class) + factor(FuelType_2001) + Distance_FS_Road + MFRI +
                              Distance_US_Highway + Distance_WUI + ERC + WindSpeed + SUP_LINE + SUP_LINE_INT +
                              LAT + LAT_LINE_INT + DIST_LAT + PREV_BURN_10Y |                    ## Other controls
                              direction_fire_FE + distance_fire_FE,                             ## FEs
                            cluster = ~FIRE_ID,
                            weights = ~Grid_Acres,
                            data = data.matched)

## What are the number of matches?

length(mgens$index.treated)  # Number of treated units matched
length(mgens$index.control)  # Number of control units matched

## Estimate using did_imputation - i.e. borusyak et al. (2024) estimator

es <- did_imputation(
  data = data.matched, yname = "BURN", gname = "L_TAU",
  tname = "distance_bin", idname = "direction_fire_FE",
  cluster_var = "FIRE_ID",
  first_stage = as.formula(paste("~", first_stage_formula)),  # Use as.formula
  wname = "Grid_Acres",
  # event-study
  horizon = TRUE,
  pretrends = -4:-1
)

es_filt <- es %>%
  filter(term %in% seq(-4, 4, 1)) %>%
  mutate(term = as.character(term)) %>%  # Convert term to character for compatibility
  bind_rows(tibble(lhs = "BURN", term = "-5", estimate = 0, std.error = 0, conf.low = 0, conf.high = 0)) %>%
  mutate(term = as.numeric(term))  # Convert to numeric for correct ordering

# Calculate the number of observations
N_obs_bin <- data.matched %>%
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

plot_matching_burn <- ggplot(es_filt, aes(x = term, y = estimate*100)) +
  # Main points
  geom_point(size = 3, color = "#0072B2") +
  # Error bars
  geom_errorbar(aes(ymin = conf.low*100, ymax = conf.high*100),
                width = 0.15, linewidth = 0.8, color = "#0072B2") +
  # Dashed reference line at 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  # Dashed vertical line at -0.5
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "black", linewidth = 0.8) +
  # Observation count above each point
  # geom_text(aes(y = conf.high*100 + 0.5, label = N_obs),
  #           size = 4, color = "black", fontface = "bold") +

  # Axis styling
  scale_x_continuous(breaks = seq(min(es_filt$term), max(es_filt$term), 1),
                     labels = interval_labels,
                     name = "Distance from treatment interaction (km)") +
  scale_y_continuous(name = "Treatment Effect (% points)") +
  labs(title = "Matching Spatial DiD - Probability of Fire Spread") +

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
  )


balance_plot_mod <- balance_plot +
  theme(plot.margin = margin(3, 14, 3, 3))   # top, right, bottom, left

plot_matching_burn_mod <- plot_matching_burn +
  theme(plot.margin = margin(3, 3, 3, 14))

## layout of Figure S9

final_layout <- (balance_plot_mod | plot_matching_burn_mod) +
  plot_layout(widths = c(1, 1.5)) +
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

## Save Figure S9

ggsave(here("output", "figures", "FigureS9.pdf"), final_layout,
       width  = 7.24,          # Science 3-column width
       height = 6,    # ~4.5 in
       units  = "in",
       dpi    = 300)

#### Save matching coefficient estimates for later tables

new_coefs <- es_filt %>%
  filter(term %in% 0:4) %>%
  arrange(term)

coef_names <- paste0("distance_bin::", new_coefs$term)

new_coef_table <- coeftable(event_twfe_Matching)

new_coef_table[coef_names, "Estimate"] <- new_coefs$estimate
new_coef_table[coef_names, "Std. Error"] <- new_coefs$std.error
new_coef_table[coef_names, "Pr(>|t|)"] <- 2 * pt(abs(new_coefs$estimate / new_coefs$std.error),
                                                 df = df.residual(event_twfe_Matching),
                                                 lower.tail = FALSE)

event_twfe_Matching$coeftable <- new_coef_table



######### Robustness Check 1 - Using Different Directions  - Creates Table S3 ###########

# Function to process the Grids_df dataset
process_grids_df <- function(file_name) {
  Grids_df <- read_csv(here("data", "intermediate", file_name))

  Grids_df <- Grids_df %>%
    mutate(
      dist_treated = ifelse(Treated_Dir == 0, 1000, L_TAU) # For Sun & Abraham Estimator
    )

  # Define treatment and control observations: post-treatment observations are considered treated if they are inside of a treatment (Treated) or directly after one (Treated_LAG)
  Grids_df_treated <- Grids_df %>%
    filter(time_to_treat > 0 & (Treated == 1 | Treated_LAG == 1))

  Grids_df_yet_treated <- Grids_df %>%
    filter(time_to_treat <= 0)

  # Combine the filtered data
  Grids_df <- bind_rows(Grids_df_treated, Grids_df_yet_treated) %>%
    filter(BURN_LAG == 1)

  # Final filtering step
  Grids_df_filt <- Grids_df %>%
    filter(time_to_treat %in% seq(-5, 4, 1) & Burn_Right_Away == 0)

  return(Grids_df_filt)
}

# Apply function to different datasets - treatment thresholds change slightly to account for the change in size of grids when lowering the number of directions.
Grids_df_filt_1_L36 <- process_grids_df("SpatialDiD_Grids_L36_K05.csv")
Grids_df_filt_1_L24 <- process_grids_df("SpatialDiD_Grids_L24_K05.csv")
Grids_df_filt_1_L18 <- process_grids_df("SpatialDiD_Grids_L18_K05_40p.csv")
Grids_df_filt_1_L12  <- process_grids_df("SpatialDiD_Grids_L12_K05_30p.csv")

###### First estimate using Fixest for saving compatability purposes with etable

run_fixest <- function(data) {
  feols(BURN ~ sunab(dist_treated, distance_bin) +
          Slope + Elev + TRI + factor(Aspect_Class) + factor(FuelType_2001) + Distance_FS_Road + MFRI +
          Distance_US_Highway + Distance_WUI + ERC + WindSpeed + FM1000 + PREV_BURN_10Y + LAT + LAT_LINE_INT +
          LOG_DIST_LAT + WIND_DIFF + DELTA_MTT + NA_DELTA_MTT + LOG_FIRE_INTENSITY + ROAD + WUI + Wilderness |
          direction_fire_FE + distance_fire_FE,
        cluster = ~FIRE_ID,
        weights = ~Grid_Acres,
        data = data)
}

# Run for each dataset
datasets <- list(Grids_df_filt_1_L36, Grids_df_filt_1_L24, Grids_df_filt_1_L18, Grids_df_filt_1_L12)
event_twfe <- lapply(datasets, run_fixest)
names(event_twfe) <- paste0("event_twfe_", 1:4)

###### Now estimate using Borusyak method

run_did_imputation <- function(data) {
  es <- did_imputation(
    data = data, yname = "BURN", gname = "L_TAU",
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

# Run for each dataset
es_results <- lapply(datasets, run_did_imputation)

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

# Apply the update function to each fixest model
event_twfe_updated <- mapply(update_fixest_coefs, event_twfe, es_results, SIMPLIFY = FALSE)
names(event_twfe_updated) <- paste0("event_twfe_", 1:4)


event_twfe_1_up <- event_twfe_updated[[1]]
event_twfe_2_up <- event_twfe_updated[[2]]
event_twfe_3_up <- event_twfe_updated[[3]]
event_twfe_4_up <- event_twfe_updated[[4]]

var_dict <- c(
  BURN = "Probability of Fire Spread",
  BURN_SEV0 = "Conditional Burn Severity",
  direction_fire_FE = "Direction-Fire",
  distance_fire_FE = "Distance-Fire",
  FIRE_ID = "Fire",
  distance_bin = "Distance Bin",
  FIRE_DAY = "Fire-Day",
  `distance_bin::0` = "$Treat_0$",
  `distance_bin::1` = "$Treat_1$",
  `distance_bin::2` = "$Treat_2$",
  `distance_bin::3` = "$Treat_3$",
  `distance_bin::4` = "$Treat_4$",
  `distance_bin::5` = "$Treat_5$"
)

setFixest_dict(var_dict)

spatDiDmods <- list()
spatDiDmods[["Burn Probability - (1)"]] <- event_twfe_1_up
spatDiDmods[["Burn Probability - (2)"]] <- event_twfe_2_up
spatDiDmods[["Burn Probability - (3)"]] <- event_twfe_3_up
spatDiDmods[["Burn Probability - (4)"]] <- event_twfe_4_up

robust_results <- etable(spatDiDmods,
                         keep = c("%distance_bin::0", "%distance_bin::1", "%distance_bin::2", "%distance_bin::3", "%distance_bin::4"),  # Use original names with "%"
                         digits = "r3",
                         digits.stats = 2,
                         fitstat = c("n", "r2"),
                         style.tex = style.tex("aer",
                                               fixef.suffix = " FEs",  # Remove the " FEs" suffix
                                               fixef.where = "var", # Keep fixed effects where they are placed
                                               yesNo = c("Yes", "No")
                         ),
                         extralines = list(
                           "No. Directions" = c("36", "24", "18", "12")
                         ),
                         signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.1), # Custom stars
                         tex = TRUE)
robust_results

robust_results %>% write_lines(here("output", "tables", "TableS3.tex"))

######### Robustness Check 2 - Testing Parallel Trends Assumption  - Creates Table S2 ###########

## Will show estimates from: i) Baseline, ii) Matched Sample, iii) No Adjacent Directions, iv) Yet-to-be treated variation

event_twfe_baseline <- event_twfe_2_up # already estimated our baseline
event_twfe_Matching <- event_twfe_Matching # already estimated our matching sample

#### Remove Adjacent Untreated Directions

## For each fire get only the adjacent "Control" directions

Grids_df <- read_csv(here("data", "intermediate", "SpatialDiD_Grids_L24_K05.csv"))

Grids_df <- Grids_df %>%
  mutate(
    Treat_Post = ifelse(Treated_Dir == 1 & time_to_treat >= 0, 1, 0),
    dist_treated = ifelse(Treated_Dir == 0, 1000, L_TAU) # For Sun & Abraham Estimator
  )

fires <- unique(Grids_df$FIRE_ID)
n_directions <- 24
Grids_df_no_adj_dir <- Grids_df[0, ]

for (fire in fires){

  Grids_df_fire <- filter(Grids_df, FIRE_ID == fire)
  treated_dirs <- unique(filter(Grids_df_fire, Treated_Dir == 1)$direction)
  control_dirs <- unique(filter(Grids_df_fire, Treated_Dir == 0)$direction)

  # Create the set of all adjacent directions with wrap-around logic
  adjacent_dirs_p1m1 <- unique(c(
    ifelse(treated_dirs - 1 == 0, n_directions, treated_dirs - 1),  # Handle wrap-around for 1
    ifelse(treated_dirs + 1 > n_directions, 1, treated_dirs + 1)    # Handle wrap-around for 48
  ))

  # # Create the set of all adjacent directions with wrap-around logic
  # adjacent_dirs_p2m2 <- unique(c(
  #   ifelse(treated_dirs - 2 <= 0, n_directions, treated_dirs - 2),  # Handle wrap-around for 1
  #   ifelse(treated_dirs + 2 > n_directions, 1, treated_dirs + 2)    # Handle wrap-around for 48
  # ))
  #
  # View the result
  adjacent_dirs <- unique(c(adjacent_dirs_p1m1))

  adjacent_control_dirs <- adjacent_dirs[adjacent_dirs %in% control_dirs]
  non_adjacent_control_dirs <- control_dirs[!(control_dirs %in% adjacent_dirs)]

  Grids_df_fire_filt_controls <- filter(Grids_df_fire, direction %in% non_adjacent_control_dirs)
  Grids_df_fire_filt_treat <- filter(Grids_df_fire, direction %in% treated_dirs)

  Grids_df_fire_filt <- rbind(Grids_df_fire_filt_controls, Grids_df_fire_filt_treat)

  Grids_df_no_adj_dir <- rbind(Grids_df_no_adj_dir, Grids_df_fire_filt)

}

Grids_df_treated <- filter(Grids_df_no_adj_dir, time_to_treat > 0) %>% filter(Treated == 1 | Treated_LAG == 1) # For any observations post-treatment make sure they are still inside of a treatment or directly after one
Grids_df_yet_treated <- filter(Grids_df_no_adj_dir, time_to_treat <= 0)

Grids_df_no_adj_dir_filt <- rbind(Grids_df_treated, Grids_df_yet_treated)
Grids_df_no_adj_dir_filt_1 <- filter(Grids_df_no_adj_dir_filt, time_to_treat %in% seq(-5, 4,1) & Burn_Right_Away == 0) %>% filter(BURN_LAG == 1)

fixest_no_adj <- run_fixest(Grids_df_no_adj_dir_filt_1)
es_no_adj <- run_did_imputation(Grids_df_no_adj_dir_filt_1)

event_twfe_noadj <- update_fixest_coefs(fixest_no_adj, es_no_adj)

#### yet-to-be treated variation

Grids_df_filt_1_L24_treated_dirs <- filter(Grids_df_filt_1_L24, Treated_Dir == 1)

fixest_treated_only <- run_fixest(Grids_df_filt_1_L24_treated_dirs)

ct <- fixest_treated_only$coeftable

blank_row <- setNames(rep(NA, ncol(ct)), colnames(ct))

ct <- rbind(
  ct,
  "distance_bin::4" = blank_row
)

fixest_treated_only$coeftable <- ct

es_treated_only <- run_did_imputation(Grids_df_filt_1_L24_treated_dirs)

event_twfe_treated_only <- update_fixest_coefs(fixest_treated_only, es_treated_only)


#### Save results to table

spatDiDmods[["Burn Probability - (1)"]] <- event_twfe_baseline
spatDiDmods[["Burn Probability - (2)"]] <- event_twfe_Matching
spatDiDmods[["Burn Probability - (3)"]] <- event_twfe_treated_only
spatDiDmods[["Burn Probability - (4)"]] <- event_twfe_noadj


#### Create Table S2

robust_variation <- etable(spatDiDmods,
                           keep = c("%distance_bin::0", "%distance_bin::1", "%distance_bin::2", "%distance_bin::3", "%distance_bin::4"),  # Use original names with "%"
                           digits = "r3",
                           digits.stats = 2,
                           fitstat = c("n", "r2"),
                           style.tex = style.tex("aer",
                                                 fixef.suffix = " FEs",  # Remove the " FEs" suffix
                                                 fixef.where = "var", # Keep fixed effects where they are placed
                                                 yesNo = c("Yes", "No")
                           ),
                           extralines = list(
                             "Baseline" = c("Yes", "No", "No", "No"),
                             "Matched Sample" = c("No", "Yes", "No", "No"),
                             "Treated Directions Only" = c("No", "No", "Yes", "No"),
                             "No Adjacent Directions" = c("No", "No", "No", "Yes")
                           ),
                           signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.1), # Custom stars
                           tex = TRUE)
robust_variation

robust_variation %>% write_lines(here("output", "tables", "TableS2.tex"))


#########   Robustness Check 3 - Alternative Estimators - Create Table S6  ###########

## i) Baseline, ii) Sun & Abraham, iii) Calloway & Santana, iv) Standard TWFEs

event_twfe_sunab <- run_fixest(Grids_df_filt_1_L24)

calloway_did <- function(data) {
  # Define the first treatment period for each grid cell
  data <- data %>%
    mutate(
      direction_fire_FE = as.numeric(as.factor(direction_fire_FE)),
      L_TAU = ifelse(is.na(L_TAU) == T, 0, L_TAU)
    )

  # Callaway & Sant'Anna estimator
  att_results <- att_gt(
    yname = "BURN",  # Outcome variable
    tname = "distance_bin",  # Time variable
    idname = "direction_fire_FE",  # Unique grid identifier
    gname = "L_TAU",  # First treatment period
    clustervars = "FIRE_ID",  # Cluster by fire event
    weightsname = "Grid_Acres",  # Weights
    allow_unbalanced_panel = T,
    xformla = ~1,
    data = data
  )

  return(att_results)
}

event_twfe_calloway <- calloway_did(Grids_df_filt_1_L24)
# summary(event_twfe_calloway)

event_twfe_calloway_dyn <- aggte(event_twfe_calloway, type = "dynamic", na.rm = T)
event_twfe_calloway_coeff <- summary(event_twfe_calloway_dyn)

calloway_estimates <- event_twfe_calloway_dyn$att.egt[5:9]
calloway_ses <- event_twfe_calloway_dyn$se.egt[5:9]


coef_names <- paste0("distance_bin::", seq(0,4,1))

new_coef_table <- coeftable(event_twfe_baseline)

new_coef_table[coef_names, "Estimate"] <- calloway_estimates
new_coef_table[coef_names, "Std. Error"] <- calloway_ses
new_coef_table[coef_names, "Pr(>|t|)"] <- 2 * pt(abs(calloway_estimates/ calloway_ses),
                                                 df = df.residual(event_twfe_baseline),
                                                 lower.tail = FALSE)

event_twfe_calloway <- event_twfe_baseline

event_twfe_calloway$coeftable <- new_coef_table

#### Standard TWFEs

twfe_standard <- function(data) {
  feols(BURN ~ i(time_to_treat, Treated_Dir, ref = -1) +
          Slope + Elev + TRI + factor(Aspect_Class) + factor(FuelType_2001) + Distance_FS_Road + MFRI +
          Distance_US_Highway + Distance_WUI + ERC + WindSpeed + FM1000 + PREV_BURN_10Y + LAT + LAT_LINE_INT +
          LOG_DIST_LAT + WIND_DIFF + DELTA_MTT + NA_DELTA_MTT + LOG_FIRE_INTENSITY + ROAD + WUI + Wilderness |
          direction_fire_FE + distance_bin,
        cluster = ~FIRE_ID,
        weights = ~Grid_Acres,
        data = data)
}

event_twfe_standard <- twfe_standard(Grids_df_filt_1_L24)
summary(event_twfe_standard)

## Change coefficient estimate names

rownames(event_twfe_standard$coeftable) <- gsub(":Treated_Dir", "", rownames(event_twfe_standard$coeftable))
rownames(event_twfe_standard$coeftable) <- gsub("time_to_treat::", "distance_bin::", rownames(event_twfe_standard$coeftable))


spatDiDmods <- list()
spatDiDmods[["Burn Probability - (1)"]] <- event_twfe_baseline
spatDiDmods[["Burn Probability - (2)"]] <- event_twfe_sunab
spatDiDmods[["Burn Probability - (3)"]] <- event_twfe_calloway
spatDiDmods[["Burn Probability - (4)"]] <- event_twfe_standard


#### Save Table S6

robust_estimators <- etable(spatDiDmods,
                            keep = c("%distance_bin::0", "%distance_bin::1", "%distance_bin::2", "%distance_bin::3", "%distance_bin::4"),  # Use original names with "%"
                            digits = "r3",
                            digits.stats = 2,
                            fitstat = c("n", "r2"),
                            style.tex = style.tex("aer",
                                                  fixef.suffix = "",  # Remove the " FEs" suffix
                                                  yesNo = c("Yes", "No")
                            ),
                            extralines = list(
                              "Baseline" = c("Yes", "No", "No", "No"),
                              "Sun & Abraham" = c("No", "Yes", "No", "No"),
                              "Callaway & Sant'anna" = c("No", "No", "Yes", "No"),
                              "Standard TWFE" = c("No", "No", "No", "Yes")
                            ),
                            signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.1), # Custom stars
                            tex = TRUE)
robust_estimators

robust_estimators %>% write_lines(here("output", "tables", "TableS6.tex"))


#########   Robustness Check 4 - Alternative Samples - Create Table S5   ###########

## Load Karen Short Fires

short <- st_read(here("data", "raw", "Short", "FPA_FOD_20221014.gdb"), layer = "Fires") %>%
  st_transform(mtbs, crs = 5070) %>%
  as.data.frame() %>%
  dplyr::select(NWCG_CAUSE_CLASSIFICATION, MTBS_ID)


mtbs <- st_read(here("data", "raw", "MTBS", "MTBS_Burn_Area", "S_USA.MTBS_BURN_AREA_BOUNDARY.shp")) %>%
  st_transform(mtbs, crs = 5070) %>%
  st_make_valid() %>%
  dplyr::rename(YEAR_MTBS = YEAR, MONTH_MTBS = STARTMONTH)

mtbs <- filter(mtbs, FIRE_TYPE == "Unknown" | FIRE_TYPE == "Wildfire") # No Wildland Fire Use Fires


mtbs_sample <- filter(mtbs, FIRE_ID %in% unique(Grids_df_filt_1_L24$FIRE_ID))

mtbs_sample <- merge(mtbs_sample, short, by.x = "FIRE_ID", by.y = "MTBS_ID", all.x = T)

mtbs_sample_filt <- filter(mtbs_sample, NWCG_CAUSE_CLASSIFICATION == "Natural")

lightning_fires <- unique(mtbs_sample_filt$FIRE_ID)

## i) Baseline, ii) Reported Ignitions, iii) USFS NF Only, iv) Remove Already Extinguished, v) Lightning Ignitions?

event_twfe_baseline <- event_twfe_2_up # already estimated our baseline

Grids_df_filt_1_L24_V2 <- process_grids_df("SpatialDiD_Grids_L24_K05_V2.csv")

Grids_df_filt_1_L24_USFS <- filter(Grids_df_filt_1_L24, USFS_NF == 1 & Wilderness == 0)

Grids_df_filt_1_L24_light <- filter(Grids_df_filt_1_L24, FIRE_ID %in% lightning_fires)

#### Create variable "Already_Extinguished", which equals one if a previous plot was unburned in a given direction

Grids_df_filt_1_L24 <- Grids_df_filt_1_L24 %>%
  arrange(direction_fire_FE, distance_bin) %>%
  group_by(direction_fire_FE) %>%
  mutate(Already_Extinguished = as.integer(lag(cummin(BURN), default = 1) == 0))

Grids_df_filt_1_L24_AE <- filter(Grids_df_filt_1_L24, Already_Extinguished == 0)



# Run for each dataset
datasets <- list(Grids_df_filt_1_L24, Grids_df_filt_1_L24_V2, Grids_df_filt_1_L24_USFS, Grids_df_filt_1_L24_AE, Grids_df_filt_1_L24_light)
event_twfe <- lapply(datasets, run_fixest)
names(event_twfe) <- paste0("event_twfe_", 1:5)

###### Now estimate using Borusyak method

datasets <- list(Grids_df_filt_1_L24, Grids_df_filt_1_L24_V2, Grids_df_filt_1_L24_AE, Grids_df_filt_1_L24_light)


# Run for each dataset
es_results <- lapply(datasets, run_did_imputation)


#### Run seperately for non-wilderness areas sample - need to get rid of Wilderness control in order to run

run_did_imputation_wilderness <- function(data) {

  controls_alt <- c("Slope", "Elev", "TRI", "Distance_FS_Road", "USFS_NF", "MFRI",
                    "Distance_US_Highway", "Distance_WUI", "ERC", "WindSpeed", "FM1000", "PREV_BURN_10Y",
                    "LAT", "LAT_LINE_INT", "LOG_DIST_LAT", "WIND_DIFF", "DELTA_MTT", "NA_DELTA_MTT",
                    "LOG_FIRE_INTENSITY", "ROAD", "WUI")
  first_stage_formula_alt <- paste(
    paste(controls_alt, collapse = " + "),  # Controls
    "|",
    paste(FEs, collapse = " + ")  # Fixed Effects
  )

  es <- did_imputation(
    data = data, yname = "BURN", gname = "L_TAU",
    tname = "distance_bin", idname = "direction_fire_FE",
    cluster_var = "FIRE_ID",
    first_stage = as.formula(paste("~", first_stage_formula_alt)),
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


es_res_1 <- run_did_imputation_wilderness(Grids_df_filt_1_L24_USFS)

es_results <- append(es_results, list(es_res_1), after = 2)

# Apply the update function to each fixest model
event_twfe_updated <- mapply(update_fixest_coefs, event_twfe, es_results, SIMPLIFY = FALSE)
names(event_twfe_updated) <- paste0("event_twfe_", 1:4)

event_twfe_baseline <- event_twfe_updated[[1]]
event_twfe_V2 <- event_twfe_updated[[2]]
event_twfe_USFS <- event_twfe_updated[[3]]
event_twfe_AE <- event_twfe_updated[[4]]
event_twfe_light <- event_twfe_updated[[5]]

spatDiDmods <- list()
spatDiDmods[["Burn Probability - (1)"]] <- event_twfe_baseline
spatDiDmods[["Burn Probability - (2)"]] <- event_twfe_V2
spatDiDmods[["Burn Probability - (3)"]] <- event_twfe_USFS
spatDiDmods[["Burn Probability - (4)"]] <- event_twfe_AE
spatDiDmods[["Burn Probability - (5)"]] <- event_twfe_light

#### Save Table S5

robust_samples <- etable(spatDiDmods,
                         keep = c("%distance_bin::0", "%distance_bin::1", "%distance_bin::2", "%distance_bin::3", "%distance_bin::4"),  # Use original names with "%"
                         digits = "r3",
                         digits.stats = 2,
                         fitstat = c("n", "r2"),
                         style.tex = style.tex("aer",
                                               fixef.suffix = " FEs",  # Remove the " FEs" suffix
                                               fixef.where = "var", # Keep fixed effects where they are placed
                                               yesNo = c("Yes", "No")
                         ),
                         extralines = list(
                           "Baseline" = c("Yes", "No", "No", "No", "No"),
                           "Reported Ignitions" = c("No", "Yes", "No", "No", "No"),
                           "USFS - Non-Wilderness" = c("No", "No", "Yes", "No", "No"),
                           "No Already Extinguished" = c("No", "No", "No", "Yes", "No"),
                           "Lightning Only" = c("No", "No", "No", "No", "Yes")
                         ),
                         signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.1), # Custom stars
                         tex = TRUE)
robust_samples

robust_samples %>% write_lines(here("output", "tables", "TableS5.tex"))


#########   Robustness Check 5 - Different Treatment Thresholds  - Create Table S7  ###########

## i) 100%, ii) 75%, iii) 50%, iv) 25%, v) Any

Grids_df_filt_1_L24_100p <- process_grids_df("SpatialDiD_Grids_L24_K05_100p.csv")
Grids_df_filt_1_L24_75p <- process_grids_df("SpatialDiD_Grids_L24_K05_75p.csv")
Grids_df_filt_1_L24_25p <- process_grids_df("SpatialDiD_Grids_L24_K05_25p.csv")
Grids_df_filt_1_L24_zero  <- process_grids_df("SpatialDiD_Grids_L24_K05_zero.csv")


## Start w/fixest to save coefficients

# Run for each dataset
datasets <- list(Grids_df_filt_1_L24_100p, Grids_df_filt_1_L24_75p,
                 Grids_df_filt_1_L24_25p, Grids_df_filt_1_L24_zero)

event_twfe <- lapply(datasets, run_fixest)
names(event_twfe) <- paste0("event_twfe_", 1:4)


## Now estimate using Borusyak method

# Run for each dataset
es_results <- lapply(datasets, run_did_imputation)

# Apply the update function to each fixest model
event_twfe_updated <- mapply(update_fixest_coefs, event_twfe, es_results, SIMPLIFY = FALSE)
names(event_twfe_updated) <- paste0("event_twfe_", 1:4)


event_twfe_1_up <- event_twfe_updated[[1]]
event_twfe_2_up <- event_twfe_updated[[2]]
event_twfe_3_up <- event_twfe_baseline # already estimated our baseline
event_twfe_4_up <- event_twfe_updated[[3]]
event_twfe_5_up <- event_twfe_updated[[4]]

spatDiDmods <- list()
spatDiDmods[["Burn Probability - (1)"]] <- event_twfe_1_up
spatDiDmods[["Burn Probability - (2)"]] <- event_twfe_2_up
spatDiDmods[["Burn Probability - (3)"]] <- event_twfe_3_up
spatDiDmods[["Burn Probability - (4)"]] <- event_twfe_4_up
spatDiDmods[["Burn Probability - (5)"]] <- event_twfe_5_up


#### Save Table S7

robust_thresholds <- etable(spatDiDmods,
                            keep = c("%distance_bin::0", "%distance_bin::1", "%distance_bin::2", "%distance_bin::3", "%distance_bin::4"),  # Use original names with "%"
                            digits = "r3",
                            digits.stats = 2,
                            fitstat = c("n", "r2"),
                            style.tex = style.tex("aer",
                                                  fixef.suffix = " FEs",  # Remove the " FEs" suffix
                                                  fixef.where = "var", # Keep fixed effects where they are placed
                                                  yesNo = c("Yes", "No")
                            ),
                            extralines = list(
                              "Baseline" = c("No", "No", "Yes", "No", "No"),
                              "% Treated Threshold" = c("100%", "75%", "50%", "25%", ">0%")
                            ),
                            signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.1), # Custom stars
                            tex = TRUE)
robust_thresholds

robust_thresholds %>% write_lines(here("output", "tables", "TableS7.tex"))



#########   Robustness 6 - Impact of Suppression Controls  - Create Table S4 ###########

## i) Baseline w/LAT, ii) Baseline no LAT, iii) Filtered Sample No Sup Lines & No Lat, iv) Filtered Sample w/LAT, v) Filtered Sample w/LAT & Lines, vi) Filter Sample Controls & Treatment only effort == 1

Grids_df_filt_1_L24_lines <- filter(Grids_df_filt_1_L24, SUP_LINES_PRESENT == 1)

Grids_df_filt_1_L24_lines <- mutate(Grids_df_filt_1_L24_lines,  SUP_LINE_TREAT = ifelse(Treated_Dir == 1 & SUP_LINE == 1, 1, 0),
                                    EFFORT = ifelse(DIST_LAT <= .20 | SUP_LINE_DIST <= .20, 1, 0))

Grids_df_filt_1_L24_lines_effort <- filter(Grids_df_filt_1_L24_lines, EFFORT == 1)


N_Fires_Baseline <- 285
N_Fires_Lines <- length(unique(Grids_df_filt_1_L24_lines$FIRE_ID))
N_Fires_Lines_Effort <- length(unique(Grids_df_filt_1_L24_lines_effort$FIRE_ID))


# Run for each dataset
datasets <- list(Grids_df_filt_1_L24, Grids_df_filt_1_L24, Grids_df_filt_1_L24_lines,
                 Grids_df_filt_1_L24_lines, Grids_df_filt_1_L24_lines, Grids_df_filt_1_L24_lines_effort)
event_twfe <- lapply(datasets, run_fixest)
names(event_twfe) <- paste0("event_twfe_", 1:6)

run_did_imputation_no_lat <- function(data) {

  controls_alt <- c("Slope", "Elev", "TRI", "Distance_FS_Road", "USFS_NF", "Wilderness", "MFRI",
                    "Distance_US_Highway", "Distance_WUI", "ERC", "WindSpeed", "FM1000", "PREV_BURN_10Y",
                    "WIND_DIFF", "DELTA_MTT", "NA_DELTA_MTT",
                    "LOG_FIRE_INTENSITY", "ROAD", "WUI")

  first_stage_formula_alt <- paste(
    paste(controls_alt, collapse = " + "),  # Controls
    "|",
    paste(FEs, collapse = " + ")  # Fixed Effects
  )


  es <- did_imputation(
    data = data, yname = "BURN", gname = "L_TAU",
    tname = "distance_bin", idname = "direction_fire_FE",
    cluster_var = "FIRE_ID",
    first_stage = as.formula(paste("~", first_stage_formula_alt)),
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

run_did_imputation_full_controls <- function(data) {

  controls_full <- c("Slope", "Elev", "TRI", "Distance_FS_Road", "USFS_NF", "Wilderness", "MFRI",
                     "Distance_US_Highway", "Distance_WUI", "ERC", "WindSpeed", "FM1000", "PREV_BURN_10Y",
                     "LAT", "LAT_LINE_INT", "LOG_DIST_LAT", "WIND_DIFF", "DELTA_MTT", "NA_DELTA_MTT",
                     "LOG_FIRE_INTENSITY", "ROAD", "WUI",  "SUP_LINE", "SUP_LINE_INT", "LOG_SUP_LINE_DIST")

  first_stage_formula_alt <- paste(
    paste(controls_full, collapse = " + "),  # Controls
    "|",
    paste(FEs, collapse = " + ")  # Fixed Effects
  )


  es <- did_imputation(
    data = data, yname = "BURN", gname = "L_TAU",
    tname = "distance_bin", idname = "direction_fire_FE",
    cluster_var = "FIRE_ID",
    first_stage = as.formula(paste("~", first_stage_formula_alt)),
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

es_res_baseline <- run_did_imputation(Grids_df_filt_1_L24)
es_res_no_lat <- run_did_imputation_no_lat(Grids_df_filt_1_L24)
es_res_lines_nocontrols <- run_did_imputation_no_lat(Grids_df_filt_1_L24_lines)
es_res_lines_lat <- run_did_imputation(Grids_df_filt_1_L24_lines)
es_res_lines_lat_lines <- run_did_imputation_full_controls(Grids_df_filt_1_L24_lines)
es_res_effort_only <- run_did_imputation_full_controls(Grids_df_filt_1_L24_lines_effort)


es_results <- list(es_res_baseline, es_res_no_lat, es_res_lines_nocontrols,
                   es_res_lines_lat, es_res_lines_lat_lines, es_res_effort_only)

# Apply the update function to each fixest model
event_twfe_updated <- mapply(update_fixest_coefs, event_twfe, es_results, SIMPLIFY = FALSE)
names(event_twfe_updated) <- paste0("event_twfe_", 1:6)


event_twfe_baseline <- event_twfe_updated[[1]]
event_twfe_no_lat <- event_twfe_updated[[2]]
event_twfe_lines_nocontrols <- event_twfe_updated[[3]]
event_twfe_lines_lat <- event_twfe_updated[[4]]
event_twfe_lines_lat_lines <- event_twfe_updated[[5]]
event_twfe_lines_effort_only <- event_twfe_updated[[6]]

spatDiDmods <- list()
spatDiDmods[["Burn Probability - (1)"]] <- event_twfe_baseline
spatDiDmods[["Burn Probability - (2)"]] <- event_twfe_no_lat
spatDiDmods[["Burn Probability - (3)"]] <- event_twfe_lines_nocontrols
spatDiDmods[["Burn Probability - (4)"]] <- event_twfe_lines_lat
spatDiDmods[["Burn Probability - (5)"]] <- event_twfe_lines_lat_lines
spatDiDmods[["Burn Probability - (6)"]] <- event_twfe_lines_effort_only

#### Save Table S4

robust_sup_controls <- etable(spatDiDmods,
                              keep = c("%distance_bin::0", "%distance_bin::1", "%distance_bin::2", "%distance_bin::3", "%distance_bin::4"),  # Use original names with "%"
                              digits = "r3",
                              digits.stats = 2,
                              fitstat = c("n", "r2"),
                              style.tex = style.tex("aer",
                                                    fixef.suffix = " FEs",  # Remove the " FEs" suffix
                                                    fixef.where = "var", # Keep fixed effects where they are placed
                                                    yesNo = c("Yes", "No")
                              ),
                              extralines = list(
                                "Baseline" = c("Yes", "No", "No", "No", "No", "No"),
                                "LAT Controls" = c("Yes", "No", "No", "Yes", "Yes", "Yes"),
                                "Suppression Line Controls" = c("No", "No", "No", "No", "Yes", "Yes"),
                                "Full Sample" = c("Yes", "Yes", "No", "No", "No", "No"),
                                "Suppression Line Sample" = c("No", "No", "Yes", "Yes", "Yes", "Yes"),
                                "Effort Only Controls" = c("No", "No", "No", "No", "No", "Yes"),
                                "No. Fires" = c(N_Fires_Baseline, N_Fires_Baseline, N_Fires_Lines,
                                                N_Fires_Lines, N_Fires_Lines_Effort, N_Fires_Lines_Effort)
                              ),
                              signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.1), # Custom stars
                              tex = TRUE)
robust_sup_controls

robust_sup_controls %>% write_lines(here("output", "tables", "TableS4.tex"))



#########   Robustness 7 - Different Event Windows - Create Table S8  ###########

process_grids_event_df <- function(file_name, event_window) {
  Grids_df <- read_csv(here("data", "intermediate", file_name))

  Grids_df <- Grids_df %>%
    mutate(
      Treat_Post = ifelse(Treated_Dir == 1 & time_to_treat >= 0, 1, 0),
      dist_treated = ifelse(Treated_Dir == 0, 1000, L_TAU) # For Sun & Abraham Estimator
    )

  # Filter treated and yet-to-be-treated observations
  Grids_df_treated <- Grids_df %>%
    filter(time_to_treat > 0 & (Treated == 1 | Treated_LAG == 1))

  Grids_df_yet_treated <- Grids_df %>%
    filter(time_to_treat <= 0)

  # Combine the filtered data
  Grids_df <- bind_rows(Grids_df_treated, Grids_df_yet_treated) %>%
    filter(BURN_LAG == 1)

  # Final filtering step
  Grids_df_filt <- Grids_df %>%
    filter(time_to_treat %in% event_window & Burn_Right_Away == 0)

  return(Grids_df_filt)
}

# Apply function to different datasets
Grids_df_filt_1_L24_25km <- process_grids_event_df("SpatialDiD_Grids_L24_K05.csv", seq(-5,4,1))
Grids_df_filt_1_L24_4km <- process_grids_event_df("SpatialDiD_Grids_L24_K05.csv", seq(-8,7,1))
Grids_df_filt_1_L24_7km <- process_grids_event_df("SpatialDiD_Grids_L24_K05.csv", seq(-15,14,1))
Grids_df_filt_1_L24_14km <- process_grids_event_df("SpatialDiD_Grids_L24_K05.csv", seq(-30,29,1))


# default results in fixest
datasets <- list(Grids_df_filt_1_L24_25km, Grids_df_filt_1_L24_4km, Grids_df_filt_1_L24_7km, Grids_df_filt_1_L24_14km)
event_twfe <- lapply(datasets, run_fixest)
names(event_twfe) <- paste0("event_twfe_", 1:4)


# Run imputation each dataset
es_results <- lapply(datasets, run_did_imputation)


# Apply the update function to each fixest model
event_twfe_updated <- mapply(update_fixest_coefs, event_twfe, es_results, SIMPLIFY = FALSE)
names(event_twfe_updated) <- paste0("event_twfe_", 1:4)

event_twfe_1_up <- event_twfe_updated[[1]]
event_twfe_2_up <- event_twfe_updated[[2]]
event_twfe_3_up <- event_twfe_updated[[3]]
event_twfe_4_up <- event_twfe_updated[[4]]

rownames(event_twfe_3_up$coeftable) <- gsub("distance_bin::(1[0-9])", "new_name_for_\\1", rownames(event_twfe_3_up$coeftable))
rownames(event_twfe_4_up$coeftable) <- gsub("distance_bin::(1[0-9])", "new_name_for_\\1", rownames(event_twfe_4_up$coeftable))
rownames(event_twfe_4_up$coeftable) <- gsub("distance_bin::([2][0-9])", "new_name_for_\\1", rownames(event_twfe_4_up$coeftable))


spatDiDmods <- list()
spatDiDmods[["Burn Probability - (1)"]] <- event_twfe_1_up
spatDiDmods[["Burn Probability - (2)"]] <- event_twfe_2_up
spatDiDmods[["Burn Probability - (3)"]] <- event_twfe_3_up
spatDiDmods[["Burn Probability - (4)"]] <- event_twfe_4_up

#### Save Table S8

robust_event_window <- etable(spatDiDmods,
                              keep = c("%distance_bin::0", "%distance_bin::1", "%distance_bin::2", "%distance_bin::3", "%distance_bin::4"),  # Use original names with "%"
                              digits = "r3",
                              digits.stats = 2,
                              fitstat = c("n", "r2"),
                              style.tex = style.tex("aer",
                                                    fixef.suffix = " FEs",  # Remove the " FEs" suffix
                                                    fixef.where = "var", # Keep fixed effects where they are placed
                                                    yesNo = c("Yes", "No")
                              ),
                              extralines = list(
                                "Baseline" = c("Yes", "No", "No", "No"),
                                "Event Window" = c("2.5 km", "4 km", "7 km", "14 km")
                              ),
                              signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.1), # Custom stars
                              tex = TRUE)
robust_event_window

robust_event_window %>% write_lines(here("output", "tables", "TableS8.tex"))

toc()

