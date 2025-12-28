# Analysis Workflow

This document describes the end-to-end workflow used to construct the analysis dataset and reproduce the empirical results in *“Wildfire damages and the cost-effective role of forest fuel treatments.”* The pipeline proceeds in sequential stages, beginning with sample definition and raw data ingestion, followed by construction of key fire progression, smoke, and simulation inputs, assembly of the Spatial DiD analysis panel, and finally estimation and robustness analyses.

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
  Intersects MTBS fire perimeters with FACTS treatment polygons to identify fires that intersect at least one completed fuel treatment prior to ignition.
  
**Key outputs:**
- List of fires included in the estimation sample

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
- Wildland Fire Emissions Inventory System (WFEIS)
- Smoke exposure estimates from Wen et al. (2023)

**Scripts:**
- `code/03_smoke/01_download_wfeis_emissions.py`  
  Downloads and extracts fire-level emissions data from WFEIS.
  
- `code/03_smoke/02_extract_smoke_exposure_wen2023.R`  
  Merges fire-level emissions with population-weighted smoke exposure estimates and imputes exposure where necessary.

**Key outputs:**
- Fire-level CO₂ and PM₂.₅ emissions
- Fire-level smoke exposure and health impact measures

---

## Step 4: FlamMap / Minimum Travel Time (MTT) Simulations

**Purpose:**  
Generate model-based predictions of fire spread behavior used to control for predictable fire dynamics in the empirical analysis.

**Key inputs:**
- LANDFIRE vegetation and fuels
- Weather conditions at ignition
- Fire ignition locations

**Scripts:**
- `code/04_mtt/01_build_mtt_landscapes.R`  
  Constructs landscape inputs required for FlamMap Minimum Travel Time (MTT) simulations.
  
- `code/04_mtt/02_run_flammap_mtt.sh`  
  Executes MTT simulations for fires in the sample.
  
- `code/04_mtt/03_process_mtt_outputs.R`  
  Processes simulation outputs to extract fire arrival times and fireline intensity measures.

**Key outputs:**
- Plot-level simulated fire arrival times
- Plot-level fireline intensity measures

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
