##########################################    main.R   ##########################################

################# Purpose: 

rm(list=ls())

# Set Path
here::i_am("code/main.R")

## Step 1 - Define sample of fires by intersecting FACTS & MTBS

source(here("code", "01_sample", "01_define_fire_sample_mtbs_facts.R"))

# ## Step 2 - Define sample of fires by intersecting FACTS & MTBS
# 
# source(here("code", "02_fire_progression", "01_define_fire_sample_mtbs_facts.R"))

## Step 3 - smoke

source(here("code", "03_smoke", "01_define_fire_sample_mtbs_facts.R"))



