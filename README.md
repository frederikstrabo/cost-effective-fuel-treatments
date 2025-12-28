# Wildfire damages and the cost-effective role of forest fuel treatments
This repository contains the code and documentation required to reproduce the main empirical analyses, tables, and figures for the paper “Wildfire damages and the cost-effective role of forest fuel treatments” (Science).

## Overview

This study integrates high-resolution spatial data on wildfires, fuel treatments, suppression effort, and economic damages across the western United States to estimate the causal effects and cost-effectiveness of forest fuel treatments. Using a quasi-experimental research design that exploits variation in the direction and distance at which fires encounter treatments, the analysis quantifies how fuel treatments reduce wildfire spread, burn severity, and downstream economics damages, including structure loss, CO₂ emissions, and PM₂.₅ exposure. The code in this repository implements the full workflow from raw data ingestion and processing through statistical estimation and figure generation. Due to size and licensing constraints the raw input datasets are not stored directly in this repository; instead, detailed instructions are provided for obtaining the required data in `/data/README.md` with the provided code reconstructing all intermediate analysis files.

A detailed description of the analysis workflow and the role of each script is provided in `code/PIPELINE.md`.
