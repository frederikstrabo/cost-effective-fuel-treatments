# Analysis Workflow

This document describes the end-to-end workflow used to construct the analysis dataset and reproduce the empirical results in *“Wildfire damages and the cost-effective role of forest fuel treatments.”* The pipeline proceeds in sequential stages, beginning with sample definition (1), followed by construction of daily fire progression (2), smoke (3), MTT simulation (4), construction of the Spatial DiD analysis panel (5), and finally the main analysis (6) and robustness (7) exercises.

Each stage produces well-defined intermediate outputs that are used as inputs to subsequent steps. Scripts are organized into numbered folders reflecting the order in which they should be executed.

# Overview of the Pipeline


## Step 1: Define Fire Sample and Treatment Intersections

**Purpose:**  
Identify the set of wildfires that intersect U.S. Forest Service fuel treatments and define the treatment–fire interactions that form the basis of the analysis sample.

**Key inputs:**
- MTBS wildfire perimeters (2017–2023)
- USFS FACTS hazardous fuel treatment polygons.

**Scripts:**
- `code/01_sample/01_define_fire_sample_mtbs_facts.R`  
  Intersects MTBS fire perimeters with FACTS treatment polygons to identify fires that intersect at least one completed fuel treatment up to 10 years prior to ignition.
  
**Key outputs:**
- List of fires included in the estimation sample - saved as `FACTS_MTBS_Fire_List.csv` in `data/intermediate` folder

These outputs define the core set of fires and treated directions used throughout the remainder of the pipeline.

---

## Step 2: Construct Daily Fire Progression (“Day of Burning”)

**Purpose:**  
Reconstruct the daily progression of each fire in the sample in order to determine when fire spread reaches different locations and intersects fuel treatments.

**Key inputs:**
- Fire perimeters
- Satellite-derived burn timing data

**Scripts:**
- `code/02_fire_progression/01_build_day_of_burning_surfaces.R`  
  Constructs interpolated daily burn surfaces for each fire in the sample.
  
**Key outputs:**
- Fire-specific day-of-burning rasters or spatial objects.
---

## Step 3: Smoke Exposure and Emissions Data

**Purpose:**  
Assemble fire-level emissions and smoke exposure measures used to quantify air quality damages.

**Key inputs:**
- Wildland Fire Emissions Inventory System (WFEIS) online calculator: [WFEIS, calculator](https://wfeis.mtri.org/calculator) 
- Smoke exposure estimates from [Wen et al. (2023), public replication code](https://github.com/jeffwen/smoke_linking_public?tab=readme-ov-file)

**Scripts:**
- `code/03_smoke/01_download_wfeis_emissions.py`  
  Downloads and extracts fire-level emissions data from WFEIS.
  
- `code/03_smoke/02_extract_smoke_exposure_wen2023.R`  
  Uses data and code from Wen et al. (2023) to get population-day weighted PM2.5 smoke exposure estimates for fires in our sample

**Key outputs:**
- Fire-level CO₂ and PM₂.₅ emissions from WFEIS - saved as `WFEIS_data.csv` in `data/raw/Smoke/WFEIS`
- Fire-level smoke exposure and health impact measures - saved as `Smoke_Fire_Effects_Wen2023.csv` saved in `data/intermediate/Wen2023`

---

## Step 4: FlamMap / Minimum Travel Time (MTT) Simulations

**Purpose:**  
Generate model-based predictions of fire spread behavior used to control for predictable fire dynamics in the empirical analysis.

**Key inputs:**
- LANDFIRE vegetation and fuels: FBFM40, CC, CH, CBH, and CBD from LANDFIRE 2001.
- Weather conditions at ignition: Wind speed and wind direction at day of ignition from gridMET.
- Fire ignition locations

**Steps:**
1. Run `code/04_mtt/01_build_mtt_landscapes.R`  
  Creates landscape files, fire ignition shapefiles, and input files for each fire in our analysis to be used in MTT. Command files are created to execute the MTT simulations in FlamMap.

2. Open FlamMap software:
  i) Take fire landscape files from `data/intermediate/MTT_Landscapes` and convert them to .lcp files and save them to `data/intermediate/FB/TestMTT/MTT_Inputs`. Rename them just replacing .tif with .lcp - i.e. save `fire1_landscape.tif` as `fire1_landscape.lcp`.
  ii) Run the simulations in FlamMap using `simulate.bat` saved in `data/intermediate/FB/TestMTT/MTT_Inputs`.
  - Outputs from simulations are shapefiles of fire arrival time and intensity:  `fireN_ArrivalTime.shp`, `fireN_INTENSITY` saved in `data/intermediate/FB/TestMTT/MTT_Inputs/MTT_Output`.
  

**Key outputs:**
- 150m simulated fire arrival times saved in `data/intermediate/FB/TestMTT/MTT_Inputs/MTT_Output`.
- 150m simulated fireline intensity saved in `data/intermediate/FB/TestMTT/MTT_Inputs/MTT_Output`.

---

## Step 5: Construct Spatial DiD Analysis Panel

**Purpose:**  
Build the plot-level panel dataset used to estimate the causal effects of fuel treatments on wildfire spread and severity.

**Key inputs:**
- Fire progression data
- Treatment intersections
- Smoke and emissions data
- MTT simulation outputs
- Weather, topography, suppression effort, and infrastructure covariates

**Scripts:**
- `code/05_panel/01_build_plot_panel.R`  
  Constructs radial plot grids for each fire, merges covariates, and assembles the final analysis-ready panel.
  
- `code/05_panel/functions_build_spatial_did.R`  
  Contains helper functions used to construct Spatial DiD plots and merge inputs.

**Key outputs:**
- Plot-level Spatial DiD panel dataset used in all estimation and robustness analyses

---

## Step 6: Main Analysis and Figure/Table Generation

**Purpose:**  
Estimate treatment effects, quantify avoided damages, and generate all figures and tables reported in the main text.

**Scripts:**
- Located in `code/06_analysis/`
- Scripts estimate treatment effects on fire spread and severity, compute cumulative impacts, perform benefit–cost calculations, and generate publication figures and tables.

---

## Step 7: Robustness and Sensitivity Analyses

**Purpose:**  
Conduct robustness checks and placebo tests reported in the Supplementary Materials.

**Scripts:**
- Located in `code/07_robustness/`
- Includes alternative specifications, placebo analyses, matching exercises, and sensitivity checks.

---

## Notes on Reproducibility

Due to the size of several datasets and external software requirements (e.g., FlamMap), some intermediate files may be generated locally and are not stored in the GitHub repository. Detailed data access instructions are provided in `data/README.md`. Users with access to all required inputs should be able to reproduce all figures and tables by executing the scripts in the order described above.
