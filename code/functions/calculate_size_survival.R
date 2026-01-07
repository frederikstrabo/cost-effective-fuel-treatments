##########################################    Calculate Size B-C Ratio    ##########################################

################# Purpose: Calculates the benefit-cost ratio based on a size bin


calculate_size_survival <- function (Size){
  
  Grids_df_filt_1 <- mutate(Grids_df_filt_1, Size_Bin = case_when(TREAT_SIZE <= q1 ~ "Small",
                                                                      TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ "Medium",
                                                                      TREAT_SIZE > q2 ~ "Large"))
  
  Grids_df_filt_1_filt_treat <- Grids_df_filt_1 %>%
    filter(Treated_Dir == 1 & Size_Bin == Size)
  
  Grids_df_filt_1_filt_untreat <- Grids_df_filt_1 %>%
    filter(Treated_Dir == 0)
  
  Grids_df_filt_1_filt <- rbind(Grids_df_filt_1_filt_treat, Grids_df_filt_1_filt_untreat)
  
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
  
  time_window <- seq(-1,9,1)
  
  Grids_df <- mutate(Grids_df, Size_Bin = case_when(TREAT_SIZE <= q1 ~ "Small",
                                                    TREAT_SIZE >= q1 & TREAT_SIZE <= q2 ~ "Medium",
                                                    TREAT_SIZE > q2 ~ "Large",))
  
  Grids_df_treat <- Grids_df %>%
    filter((Treated_Dir == 1 & Size_Bin == Size))
  
  Grids_df_untreat <- Grids_df %>%
    filter(Treated_Dir == 0)
  
  Grids_df <- rbind(Grids_df_treat, Grids_df_untreat)
  
  treated_data <- Grids_df %>% filter(Treated_Dir == 1 & Burn_Right_Away == 0) %>% filter(time_to_treat %in% time_window)
  
  treated_data$pred_without_treatment <- predict(OLS_Pred_Model, newdata = treated_data, type = "response")
  treated_data$pred_without_treatment_BS <- predict(OLS_Pred_Model_Burn, newdata = treated_data, type = "response")
  treated_data$pred_without_treatment_MED_HIGH_BS <- predict(OLS_Pred_Model_Med_High_Burn, newdata = treated_data, type = "response")
  treated_data$pred_without_treatment_HIGH_BS <- predict(OLS_Pred_Model_High_Burn, newdata = treated_data, type = "response")
  
  
  sum(is.na(treated_data$pred_without_treatment))
  sum(is.na(treated_data$pred_without_treatment_BS))
  
  ex <- dplyr::select(treated_data, direction_fire_FE, BURN, time_to_treat, dist_treated, distance_bin, pred_without_treatment) %>% arrange(direction_fire_FE)
  
  treated_data <- filter(treated_data, is.na(pred_without_treatment) == F)
  
  ## Option 1 - for predictions > 1 set to 1, < 0, set to 0
  
  # treated_data <- mutate(treated_data, pred_without_treatment = case_when(pred_without_treatment > 1 ~ 1,
  #                                                                         pred_without_treatment < 0 ~ 0,
  #                                                                         pred_without_treatment >= 0 & pred_without_treatment <= 1 ~ pred_without_treatment))


  treated_data[is.na(treated_data$pred_without_treatment_BS), "pred_without_treatment_BS"] <- 0
  treated_data[is.na(treated_data$pred_without_treatment_MED_HIGH_BS), "pred_without_treatment_MED_HIGH_BS"] <- 0
  treated_data[is.na(treated_data$pred_without_treatment_HIGH_BS), "pred_without_treatment_HIGH_BS"] <- 0
  
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
  
  survival_plot <- ggplot(survivor_plot, aes(x = time_to_treat)) +
    # Dashed vertical line at -0.5
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_line(aes(y = Mean_Survival_Treat, color = "With Treatment"), size = 1.2) +
    geom_line(aes(y = Mean_Survival_No_Treat, color = "Without Treatment"), size = 1.2) +
    labs(
      title = "Predicted Burn Probability by Distance From Treatment",
      x = "Distance from treatment interaction (0.5km)",
      y = "Predicted Burn Probability (% pt)",
      color = "Treatment Status"
    ) +
    scale_color_manual(values = c("With Treatment (Observed)" = "#0072B2", "Without Treatment (Counterfactual)" = "#D55E00")) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),  # Remove minor gridlines for a cleaner look
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add axes
      axis.line = element_blank(),  # Avoid duplicate axis lines
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 14),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  
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
  
  survivor_plot_boot$distance <- seq(0, 2.5, 0.5)
  
  
  starting_BS <- filter(survivor_plot_boot, time_to_treat == -1)$Pct_Burned_MS_HS_Treat
  
  survivor_plot_boot[1, c("Treated_Mean", "Treated_Lower", "Treated_Upper", 
                          "Untreated_Mean", "Untreated_Lower", "Untreated_Upper")] <- rep(1, 6)
  
  survivor_plot_boot[1, c("Diff_Mean", "Diff_Lower", "Diff_Upper", 
                          "Diff_BS_Mean", "Diff_BS_Lower", "Diff_BS_Upper")] <- rep(0, 6)
  
  
  survivor_plot_boot[1, c("Treated_BS_Mean", "Treated_BS_Lower", "Treated_BS_Upper", 
                          "Untreated_BS_Mean", "Untreated_BS_Lower", "Untreated_BS_Upper")] <- rep(starting_BS, 6)
  
  
  if (Size == "Small"){
    title_text <- "Small (75-600 Acres)"
  }
  
  if (Size == "Medium"){
    title_text <- "Medium (600-2400 Acres)"
  }
  
  if (Size == "Large"){
    title_text <- "Large (> 2400 Acres)"
  }
  
  survival_plot_boot <- ggplot(survivor_plot_boot, aes(x = distance)) +
    # Confidence intervals as ribbons
    geom_ribbon(aes(ymin = Treated_Lower*100, ymax = Treated_Upper*100, fill = "With Treatment (Observed)"), alpha = 0.2) +
    geom_ribbon(aes(ymin = Untreated_Lower*100, ymax = Untreated_Upper*100, fill = "Without Treatment (Counterfactual)"), alpha = 0.2) +
    # Mean survival lines
    geom_line(aes(y = Mean_Survival_Treat*100, color = "With Treatment (Observed)"), size = 1.2) +
    geom_line(aes(y = Mean_Survival_No_Treat*100, color = "Without Treatment (Counterfactual)"), size = 1.2) +
    
    # Labels and titles
    labs(
      title = title_text,
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
  
  
  survival_MH_BS_plot_boot <- ggplot(survivor_plot_boot, aes(x = distance)) +
    # Confidence intervals as ribbons
    geom_ribbon(aes(ymin = 100*Treated_BS_Lower, ymax = Treated_BS_Upper*100, fill = "With Treatment (Observed)"), alpha = 0.2) +
    geom_ribbon(aes(ymin = 100*Untreated_BS_Lower, ymax = Untreated_BS_Upper*100, fill = "Without Treatment (Counterfactual)"), alpha = 0.2) +
    # Mean survival lines
    geom_line(aes(y = Pct_Burned_MS_HS_Treat*100, color = "With Treatment (Observed)"), size = 1.2) +
    geom_line(aes(y = Pct_Burned_MS_HS_No_Treat*100, color = "Without Treatment (Counterfactual)"), size = 1.2) +
    
    # Labels and titles
    labs(
      title = title_text,
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

  
  # treated_data <- mutate(treated_data, Acres_Burned_With_Treat = survival_with_treatment*Grid_Acres,
  #                        Acres_Burned_No_Treat = survival_without_treatment*Grid_Acres,
  #                        Acres_Burned_Realized = BURN*Grid_Acres,
  #                        Acres_Burned_HS_With_Treat = survival_with_treatment*Grid_Acres*pred_with_treatment_HIGH_BS,
  #                        Acres_Burned_HS_Without_Treat = survival_without_treatment*Grid_Acres*pred_without_treatment_HIGH_BS)
  # 
  # acres_burned_savings <- sum(treated_data$Acres_Burned_No_Treat) - sum(treated_data$Acres_Burned_With_Treat) # Savings only considering 5 leads
  # acres_burned_savings_pct <- (sum(treated_data$Acres_Burned_No_Treat) - sum(treated_data$Acres_Burned_With_Treat))/sum(treated_data$Acres_Burned_Realized)
  # 
  # treated_data <- mutate(treated_data, Structure_Burned_With_Treat = survival_with_treatment*Struc_Count,
  #                        Structure_Burned_No_Treat = survival_without_treatment*Struc_Count,
  #                        Structure_Burned_Realized = BURN*Struc_Count,
  #                        Structure_Value_Lost_With_Treat = survival_with_treatment*TOT_STRUC_VAL,
  #                        Structure_Value_Lost_No_Treat = survival_without_treatment*TOT_STRUC_VAL,
  #                        Structure_Value_Lost_Realized = BURN*TOT_STRUC_VAL)
  # 
  # house_val_savings <- sum(treated_data$Structure_Value_Lost_No_Treat, na.rm = T) - sum(treated_data$Structure_Value_Lost_With_Treat, na.rm = T)
  # 
  # treated_data <- mutate(treated_data, PM25_With_Treat = Acres_Burned_With_Treat*FIRE_PM25_PER_ACRE,
  #                        PM25_No_Treat = Acres_Burned_No_Treat*FIRE_PM25_PER_ACRE,
  #                        PM25_Realized = Acres_Burned_Realized*FIRE_PM25_PER_ACRE,
  #                        CO2_With_Treat = Acres_Burned_With_Treat*FIRE_CO2_PER_ACRE,
  #                        CO2_No_Treat = Acres_Burned_No_Treat*FIRE_CO2_PER_ACRE,
  #                        CO2_Realized = Acres_Burned_Realized*FIRE_CO2_PER_ACRE)
  # 
  # carbon_savings <- (sum(treated_data$CO2_No_Treat) - sum(treated_data$CO2_With_Treat))*185 # CO2 * Carbon price
  # 
  # treated_data <- mutate(treated_data,
  #                        Deaths_With_Treat = Acres_Burned_With_Treat*Deaths_PER_ACRE,
  #                        Deaths_No_Treat = Acres_Burned_No_Treat*Deaths_PER_ACRE,
  #                        Deaths_Realized = Acres_Burned_Realized*Deaths_PER_ACRE,
  #                        DeathCost_With_Treat = Acres_Burned_With_Treat*DeathCost_PER_ACRE,
  #                        DeathCost_No_Treat = Acres_Burned_No_Treat*DeathCost_PER_ACRE,
  #                        DeathCost_Realized = Acres_Burned_Realized*DeathCost_PER_ACRE,
  #                        EarningsLoss_With_Treat = Acres_Burned_With_Treat*EarningsLoss_PER_ACRE,
  #                        EarningsLoss_No_Treat = Acres_Burned_No_Treat*EarningsLoss_PER_ACRE,
  #                        EarningsLoss_Realized = Acres_Burned_Realized*EarningsLoss_PER_ACRE)
  # Death_savings <- (sum(treated_data$DeathCost_No_Treat, na.rm = T) - sum(treated_data$DeathCost_With_Treat, na.rm = T)) # Deaths*VSL
  # Earnings_savings <- (sum(treated_data$EarningsLoss_No_Treat, na.rm = T) - sum(treated_data$EarningsLoss_With_Treat, na.rm = T)) # Earnings Estimate
  # 
  # facts_act_grouped_06_23_filt_size <- filter(facts_act_grouped_06_23_filt_1, Size_Bin == Size)
  # 
  # Total_Cost <- sum(facts_act_grouped_06_23_filt_size$TOT_COST)
  # N_Treats <- length(unique(treated_data$TREAT_ID))
  # mu_s <- (house_val_savings + carbon_savings + Death_savings + Earnings_savings)/N_Treats # estimated benefit per treatment conditional on intersection
  # 
  # facts_act_grouped_06_23_filt_size <- mutate(facts_act_grouped_06_23_filt_size, Benefits = mu_s*lambda)
  # 
  # BCR <- sum(facts_act_grouped_06_23_filt_size$Benefits)/Total_Cost
  
  return(list(
    survivor_plot = survival_plot_boot,
    survivor_BS_plot = survival_MH_BS_plot_boot
  ))
  
}