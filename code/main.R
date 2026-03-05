##########################################    main.R   ##########################################

################# Purpose: This script runs all the scripts used in R to replicate the analysis.  
#################             See "PIPELINE.md" for more information on how to run analysis needed to replicate outside of R:
#################               i.e. running the Step 3.1 in Python or the MTT simulations in FlamMap. 


################# Note: Users who wish to reproduce the main estimation results, tables, and figures without rerunning the full 
#################         preprocessing pipeline may begin at Step 6 using the cleaned, analysis-ready datasets provided in `data/intermediate/`. 
#################         This path allows replication of all main results and figures reported in the paper. 

rm(list=ls())

pacman::p_load(here)


# Set Path
here::i_am("code/main.R")

#### Set the path to your project to activate the environment

setwd("/home/ftstrabo/Documents/cost-effective-fuel-treatments") # CHANGE THIS
source("renv/activate.R")


## Step 1: Define Fire Sample and Treatment Intersections

source(here("code", "01_sample", "01_define_fire_sample_mtbs_facts.R"))

## Step 2: Construct Daily Fire Progression (“Day of Burning”)

source(here("code", "02_fire_progression", "01_dob_interpolation.R"))

source(here("code", "02_fire_progression", "02_parks_facts_list.R"))

## Step 3: Smoke Exposure and Emissions Data

## Step 3.1 run "01_download_wfeis_emissions.py" in python

source(here("code", "03_smoke", "02_extract_smoke_exposure_wen2023.R"))

## Step 4: FlamMap / Minimum Travel Time (MTT) Simulations

source(here("code", "04_mtt", "01_build_mtt_landscapes.R"))

## See "PIPELINE.md" for instructions for running the MTT simulations in FlamMap

## Step 5: Construct Spatial DiD Analysis Panel

source(here("code", "05_panel", "01_build_plot_panel.R"))
source(here("code", "05_panel", "02_build_incomplete_panel.R")) # only run if you want to run the incomplete treatment robustness check in scrip "07_robustness/03_conditional_effects_robustness.R"

################  Analysis-only replication

## Step 6: Main Analysis and Figure/Table Generation

source(here("code", "06_analysis", "01_descriptive_stats.R"))
source(here("code", "06_analysis", "02_conditional_effects.R"))
source(here("code", "06_analysis", "03_cumulative_effects.R"))

## Step 7: Robustness and Sensitivity Analyses

source(here("code", "07_robustness", "01_cumulative_effects_yet_to_be_treated.R"))
source(here("code", "07_robustness", "02_cumulative_effects_matching.R"))
source(here("code", "07_robustness", "03_conditional_effects_robustness.R"))

## Step 8: Create maps

source(here("code", "08_maps", "01_create_maps.R"))

