# Wildfire damages and the cost-effective role of forest fuel treatments
Repository supporting [Strabo, Bryan, & Reimer (2026)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5337929).

## Overview

This study integrates high-resolution spatial data on wildfires, fuel treatments, suppression effort, and economic damages across the western United States to estimate the causal effects and cost-effectiveness of forest fuel treatments. Using a quasi-experimental research design that exploits variation in the direction and distance at which fires encounter treatments, the analysis quantifies how fuel treatments reduce wildfire spread, burn severity, and downstream economics damages, including structure loss, CO₂ emissions, and PM₂.₅ exposure. 

This GitHub repository contains the full codebase implementing the workflow from raw data ingestion and processing through statistical estimation and figure generation.

- A detailed description of the analysis workflow and the role of each script is provided in `code/PIPELINE.md`. It is recommended to read this before replication.
- A detailed description of the data raw data directory is provided in `data/data_readme.md`.


## Replication Archive (Code + Data)

Due to size constraints raw input datasets are not stored directly in this GitHub repository.

The complete replication package, including all code and datasets necessary to reproduce the results and figures reported in the paper, is archived at Zenodo:

DOI: 10.5281/zenodo.XXXXXXX

## Directory Structure

The directory structure for the project is as follows.

Directory                                  | Description
-------------------------------------------|-----------------------------------------
`data`              | Data folder with subdirectories `raw`, `intermediate`, and `temp`.
`code`              | All the main R scripts used for data cleaning, data wrangling, analysis, figures, etc.
`code/functions`    | Any custom R functions
`output`            | Outputs from our R scripts, such as plots, maps, tables, estimation results etc. are saved in this directory.



## License

The code in this replication package is licensed under the MIT License.

The underlying datasets are publicly available from the original data providers cited in the manuscript and remain subject to the terms of use specified by those providers.

## Citation

If you use this code or replication archive, please cite:

Strabo, F., Bryan, C., & Reimer, M. (2026). Wildfire damages and the cost-effective role of forest fuel treatments. Science (forthcoming).

