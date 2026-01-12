# Analysis Workflow

This document describes the end-to-end workflow used to construct the analysis dataset and reproduce the empirical results in *“Wildfire damages and the cost-effective role of forest fuel treatments.”* The pipeline proceeds in sequential stages, beginning with sample definition (1), followed by construction of daily fire progressions (2), smoke impacts (3), MTT simulations (4), construction of the spatial DiD analysis panel (5), the main analysis (6), robustness (7) exercises, and creation of maps (8).

Each stage produces well-defined intermediate outputs that are used as inputs to subsequent steps. Scripts are organized into numbered folders reflecting the order in which they should be executed.

## Platform and software notes

The full preprocessing pipeline has been developed and tested on macOS and Linux. Several steps rely on fork-based parallel processing via `parallel::mclapply()` to speed up large spatial operations; because forking is not supported on Windows, these steps will not run natively on Windows without modification.

Windows users have several options:
1. Run the pipeline in a Linux environment (recommended), such as Windows Subsystem for Linux (WSL2).
2. Modify the relevant scripts to replace `mclapply()` with a cross-platform alternative (e.g., `future.apply::future_lapply()` using a multisession plan, or serial execution via `lapply()`).

In addition, the FlamMap / Minimum Travel Time (MTT) simulations (Step 4) require the FlamMap software and some manual user interaction. To facilitate reproducibility, this repository includes the processed MTT outputs used in the analysis, allowing users to reproduce all main results without rerunning the MTT simulations. Instructions for regenerating MTT outputs are provided for users with access to FlamMap.

## Restricted data and reproducibility notes

Some specifications in the analysis incorporate information on wildfire suppression effort constructed using U.S. Forest Service large airtanker (LAT) drop location data. These data were provided by the U.S. Forest Service under restricted access and are not publicly distributable, and therefore are not included in this repository. Users with approved access to the LAT data may place the files in `data/raw/LAT/`, in which case the suppression-augmented specifications will be reproduced automatically.

The LAT data are used to construct plot-level indicators of proximity to aerial suppression effort and to examine heterogeneity in treatment effects by suppression presence. The baseline treatment effect estimates are similar without the inclusion of these controls (see Table S4). 

---

## Replication paths

This repository supports two primary replication workflows:

**1. Full pipeline replication (macOS/Linux)**  
Users with access to all required software and inputs may reproduce the analysis starting from raw data by executing Steps 1–8 in sequence. This path reconstructs all intermediate datasets used in the paper.

**2. Analysis-only replication (all platforms)**  
Users who wish to reproduce the main estimation results, tables, and figures without rerunning the full preprocessing pipeline may begin at Step 6 using the cleaned, analysis-ready datasets provided in `data/intermediate/`. This path allows replication of all main results and figures reported in the paper. 

The analysis-only replication path is recommended for users who do not have access to restricted datasets (e.g., LAT drops) or external software (e.g., FlamMap).


# Overview of the Pipeline

## Step 1: Define Fire Sample and Treatment Intersections

**Purpose:**  
Identify the set of wildfires that intersect U.S. Forest Service fuel treatments and define the treatment–fire interactions that form the basis of the analysis sample.

**Key inputs:**
- MTBS wildfire perimeters (2017–2023)
- USFS FACTS hazardous fuel treatment polygons.

**Scripts:**
- `code/01_sample/01_define_fire_sample_mtbs_facts.R`  
  - Intersects MTBS fire perimeters with FACTS treatment polygons to identify fires that intersect at least one completed fuel treatment up to 10 years prior to ignition.
  - **output:**: List of fires included in the estimation sample - saved as `FACTS_MTBS_Fire_List.csv` in `data/intermediate` folder
  - Estimated run time: ~ 25 minutes.
---

## Step 2: Construct Daily Fire Progression (“Day of Burning”)

**Purpose:**  
Reconstruct the daily progression of each fire in the sample in order to determine when fire spread reaches different locations and intersects fuel treatments.

**Key inputs:**
- Fire perimeters from MTBS. 
- Satellite-derived burn timing data from FIRMS.
- List of fires included in the estimation sample - `FACTS_MTBS_Fire_List.csv`.

**Scripts:**
- `code/02_fire_progression/01_dob_interpolation.R`  
  - Constructs interpolated daily burn surfaces for each fire in the sample.
  - **output**: folder of day of burning rasters called `Jan2025_DOB_fires` saved in `data/raw/PARKS`.
  - Estimated run time: 21 hours.

- `code/02_fire_progression/02_parks_facts_list.R`  
  - Creates a list of all the daily of burning fire perimeter rasters used in our sample. This list is used as input when we build our plot panels in `01_build_plot_panel.R`.
  - **output**: `FACTS_Parks_Fire_List.csv` saved in `data/intermediate`.
  - Estimated run time: ~1 hour and 40 minutes. 

---

## Step 3: Smoke Exposure and Emissions Data

**Purpose:**  
Assemble fire-level emissions and smoke exposure measures used to quantify air quality damages.

**Key inputs:**
- Wildland Fire Emissions Inventory System (WFEIS) online calculator: [WFEIS, calculator](https://wfeis.mtri.org/calculator) 
- Smoke exposure estimates from [Wen et al. (2023), public replication code](https://github.com/jeffwen/smoke_linking_public?tab=readme-ov-file)

**Scripts:**
- `code/03_smoke/01_download_wfeis_emissions.py`  
  - Downloads and extracts fire-level emissions data from WFEIS.
  - **output**: Fire-level CO₂ and PM₂.₅ emissions from WFEIS - saved as `WFEIS_data.csv` in `data/raw/WFEIS`.
  - Estimated run time: ~3 hours.
  
- `code/03_smoke/02_extract_smoke_exposure_wen2023.R`  
  - Uses data and code from Wen et al. (2023) to get population-day weighted PM₂.₅ smoke exposure estimates for fires in our sample.
  - **output**: Fire-level smoke exposure and health damage measures - saved as `Smoke_Fire_Effects_Wen2023.csv` saved in `data/intermediate/Wen2023`.
  - Estimated run time: ~5 minutes.

---

## Step 4: FlamMap / Minimum Travel Time (MTT) Simulations

**Purpose:**  
Generate model-based predictions of fire spread behavior used to control for predictable fire dynamics in the empirical analysis.

**Key inputs:**
- LANDFIRE vegetation and fuels: Fire Behavior Fuel Model 40 (FBFM40), Canopy Cover (CC), Canopy Height (CH), Canopy Base Height (CBH), and Canopy Bulk Denisty (CBD) from LANDFIRE 2001.
- Weather conditions at ignition: Wind speed and wind direction at day of ignition from gridMET.
- MTBS fire perimeters and day of burning rasters.

**Steps:**
1. Run `code/04_mtt/01_build_mtt_landscapes.R`  
  - Creates landscape files, fire ignition shapefiles, and input files for each fire in our analysis to be used in MTT. Command files are created to execute the MTT simulations in FlamMap.
  - **Outputs**:
    - landscape.tif files for all fires saved in `data/intermediate/MTT_Landscapes`.
    - fire ignition shapefiles saved in `data/raw/FB/TestMTT/MTT_Inputs`.
    - fire .input files saved in `data/raw/FB/TestMTT/MTT_Inputs`.
    - `simulateCMD.txt` and `simulate.bat` saved in `data/raw/FB/TestMTT/MTT_Inputs`.
  - Estimated run time: ~45 mins.

2. Open FlamMap software:
    - Take fire landscape files from `data/raw/MTT_Landscapes` and convert them to .lcp files and save them to `data/raw/FB/TestMTT/MTT_Inputs`. Rename them just replacing .tif with .lcp - i.e. save `fire1_landscape.tif` as `fire1_landscape.lcp`.
    - Run the simulations in FlamMap using `simulate.bat` saved in `data/raw/FB/TestMTT/MTT_Inputs`.
    - **Outputs** from simulations are shapefiles of fire arrival time and intensity:  `fireN_ArrivalTime.shp`, `fireN_INTENSITY` saved in `data/raw/FB/TestMTT/MTT_Inputs/MTT_Output`.
    - Estimated run time: ~2 days.

---

## Step 5: Construct Spatial DiD Analysis Panel

**Purpose:**  
Build the plot-level panel dataset used to estimate the causal effects of fuel treatments on wildfire spread and severity.

**Key inputs:**
- Fire progression & perimeter data.
- Fuel treatment data.
- PM₂.₅ smoke exposure and CO₂ emissions data.
- MTT simulation outputs.
- Weather, topography, suppression effort, and infrastructure covariates.

**Scripts:**
- `code/05_panel/01_build_plot_panel.R`  
  - For each fire in our sample we construct plots which are defined by their unique direction and distance from the fires ignition point. For each plot we calculate their associated wildfire outcomes, treatment status, and relevant control covariates to be used in our spatial DiD analysis. The script develops the function, "create_SpatialDiD_plot_panel", which allows flexibility to change the number of radial directions = L, the distance of each plot = K, which will be tested in "07_robustness".
  - **Note**: This code defaults to running the baseline plots. At the end of the script the user can uncomment code used to construct alternative samples used in "07_robustness", which greatly increases the run time.
  - **output**: The baseline plots used in spatial DiD analysis `SpatialDiD_Grids_L24_K05.csv` saved in `data/temp`. Note the panel is finished in `06_analysis/01_descriptive_stats.R` where smoke exposure is imputed for fires without PM₂.₅ exposure estimates.
  - Optional **outputs**: Different direction panels: `"SpatialDiD_Grids_L36_K05.csv`, `"SpatialDiD_Grids_L18_K05_40p.csv`, and `SpatialDiD_Grids_L12_K05_30p.csv`, reported ignition points: `SpatialDiD_Grids_L24_K05_V2.csv`, and alternative treatment thresholds: `SpatialDiD_Grids_L24_K05_100p.csv`, `SpatialDiD_Grids_L24_K05_75p.csv`, `SpatialDiD_Grids_L24_K05_25p.csv`, `SpatialDiD_Grids_L24_K05_zero.csv` all of which are saved to `data/intermediate`.
  - Estimated run time: ~X hours for baseline plots and ~Y hours with the inclusion of all alternative plots used to test robustness.
- `code/05_panel/02_build_incomplete_panel.R`.
  -   Using the similar code from `05_panel/01_build_plot_panel.R`, create a panel of plots using only incomplete fuel treatment projects for a placebo test. Code is not exactly the same as `05_panel/01_build_plot_panel.R"`to avoid unnecessary computation.
  -   **output**: `SpatialDiD_Grids_24L_K05_Incomp.csv` saved in `data/intermediate`.
  -   Estimated run time: ~X min.

---

## Step 6: Main Analysis and Figure/Table Generation

**Purpose:**  
Estimate treatment effects, quantify avoided damages, and generate all figures and tables reported in the main text.

**Scripts:**
- `code/06_analysis/01_descriptive_stats.R`  
  - Create maps and time series used to create Figure 1, impute PM2.5 smoke exposure, and calculate USFS budget imputed footprint cost/acre.
  - **outputs**:
    - `SpatialDiD_Grids_L24_K05.csv` - saved in `data/intermediate`.
    - Figure S10 saved as `Emissions_Exposure_Plot.pdf` stored in `output/figures`.
    - Figure 1 saved as `Figure1.pdf` stored in `output/figures`.
   
- `code/06_analysis/02_conditional_effects.R`.
  -  Run the baseline spatial DiD regressions, create event study plots, and heterogeneity analysis.
  -  **outputs**:
    - Figure 3 saved as `Figure3.pdf` saved in `output/figures`.
    - Figure 5 saved as `Figure5.pdf` saved in `output/figures`.
    - Figure S5 saved as `HeterogeneityBurnSeverity.pdf` saved in `output/figures`.
    - Figure S7 saved as `EventWindows.pdf` saved in `output/figures`.
    - Table S1 saved as `SpatialDiDBaseline.tex` saved in `output/tables`.
  - Estimated run time: ~2 min.
      
- `code/06_analysis/03_cumulative_effects.R`.
  -   Conduct the cumulative effects (i.e. "survival analysis", Figure 4), calculate conditional treatment benefits (Table 1, A & B) , and conduct the ex-ante benefit cost ratio (BCR) analysis (Table 1 C).
  -   **outputs**:
    - Figure 4 saved as `Figure4.pdf` saved in `output/figures`.
    - Figure S1 saved as `Survival_Difference.pdf` saved in `output/figures`.
    - Figure S3  saved as `EventStudy_Survival_High_Ordinal_BS.pdf` saved in `output/figures`.
    - Tables 1 A), B) & C) saved as `PhysicalBenefits.tex`, `EconomicBenefits.tex`, and `ExAnteBCRs.tex` saved in `output/tables`.
    - Table S9 saved as `BenefitCostRatios.tex` saved in `output/tables`.
    - Figure S6 saved as `TreatmentInteractionSurvival.pdf` saved in `output/figures`.
    - Figure S4 saved as `SurvivalSizeHeterogeneity.pdf` saved in `output/figures`.
    - `BCR_Robustness.csv`  saved in `data/temp`.      
  -   Estimated run time: ~45 min.
  

---

## Step 7: Robustness and Sensitivity Analyses

**Purpose:**  
Conduct robustness checks and placebo tests reported in the Supplementary Materials.

**Scripts:**
- Located in `code/07_robustness/`
- `code/07_robustness/01_cumulative_effects_yet_to_be_treated.R`.
  -   Re-run the survival and BCR analysis using the "yet-to-be treated" control specification.
  -   **output**: `BCR_Robustness_up.csv`  saved in `data/temp`.
  -   Estimated run time: ~20 min.
- `code/07_robustness/02_cumulative_effects_matching.R`.
  -   Re-run the survival and BCR analysis using the matched control specification.
  -   **output**: Table S10 saved as `RobustnessBCR.tex` saved in `output/tables`.
  -   Estimated run time: ~22 min.

- `code/07_robustness/03_conditional_effects_robustness.R`.
  -   Run robustness checks on the conditional effects spatial DiD regressions reported in supplemental appendix.
  -   **outputs**:
    - Figure S8 saved as `IncompletePlacebo.pdf` stored in `output/figures`.
    - Figure S9 saved as `MatchingDiD.pdf` in `output/figures`.
    - Figure S9 saved as `MatchingDiD.pdf` in `output/figures`.
    -  Table S3 saved as `SpatialDiDRobustDirections.tex` in `output/tables`.
    -  Table S2 saved as `SpatialDiDRobustParallelTrends.tex` in `output/tables`.
    -  Table S2 saved as `SpatialDiDRobustEstimators.tex` in `output/tables`.
    -  Table S5 saved as `SpatialDiDSamples.tex` in `output/tables`.
    -  Table S7 saved as `SpatialDiDRobustTreatmentThresholds.tex` in `output/tables`.
    -  Table S4 saved as `SpatialDiDSupControls.tex` in `output/tables`.
    -  Table S8 saved as `SpatialDiDRobustEventWindow.tex` in `output/tables`.
  -   Estimated run time: ~1.5 hours.

---


## Step 8: Create maps (optional)

**Purpose:**  
Creates maps of fires used to create Figure 2 and Figure S2.

**Scripts:**
- Located in `code/08_maps/01_create_maps.R`
  - Creates maps of fires and plots used to create Figure 2 and Figure S2.
  - **outputs**:
    - Figure 2 saved as `Figure2.pdf` stored in `output/figures`.
    - Figure S2 saved as `Fire-FT-Examples.pdf` in `output/figures`.
  - Estimated run time: ~ 10 min.
