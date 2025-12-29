# Data Sources

## Setting Up "data/raw" directory

In order for the code to run we need to have the same folder/file names in "data/raw" for the data downloaded from the same links: 

Data folder     | Description
--------------------|------------------------------------------
`ACS`                             | American Community Survey (ACS) Census Block Group Housing Data downloaded from: [IPUMS](https://usa.ipums.org/usa/)
`s`                                 | `ACS\nhgis0003_shape`-shapefiles of 2020 Census Blocks 
 `s`                                 | `ACS\nhgis0003_csv` - csv's of ACS 2013-2017,...,2018-2022 surveys with Median Housing Values in a census block.
`CommunitiesRisk`                 | [Wildfire Risk to Communities: Spatial datasets of wildfire risk for populated areas in the United States (2nd Edition](https://www.fs.usda.gov/rds/archive/catalog/RDS-2020-0060-2)
`FACTS`                           | USFS Fuel Treatment Data downloaded from: [USDA Forest Service FSGeodata Clearinghouse - Download National Datasets](https://data.fs.usda.gov/geodata/edw/datasets.php) under "Hazardous Fuel Treatment Reduction: Polygon" (shapefile) - Rename downloaded folder to "FACTS"
`FRED`                            | Inflation adjustment data and U.S. elderly population share taken from Federal Reserve Bank of St. Louis (FRED) 
                                  | `CPIAUCSL.csv` yearly consumer index for inflation adjustment downloaded from [FRED](https://fred.stlouisfed.org/series/CPIAUCSL)
                                  | `SPPOP65UPTOZSUSA.csv` yearly population ages 65 and above in U.S. downloaded from [FRED](https://fred.stlouisfed.org/series/SPPOP65UPTOZSUSA)
`gridMET`                         | Gridded Surface Meteorological dataset (gridMET). Yearly energy Release Component (erc), 1000 hour fuel moisture (fm1000), wind speed (vs), and wind direction (th) downloaded from [gridMET](https://www.northwestknowledge.net/metdata/data/)
`Highways`                        | Shapefile of U.S. Highways from the U.S. Census Bureau downloaded from [TIGER/Line Shapefile, 2016, nation, U.S., Primary Roads National Shapefile](https://catalog.data.gov/dataset/tiger-line-shapefile-2016-nation-u-s-primary-roads-national-shapefile)
`ICS-209`                         |
`LANDFIRE`                        | Slope, aspect, elevation, mean fire return interval (MFRI), exisiting vegetation type, fire behavior fuel model 40 (FBFM40), canopy cover (CC), canopy height (CH), canopy base height (CBH), and canopy bulk denisty (CBD) for CONUS downloaded from [LANDFIRE](https://www.landfire.gov/data)
`LAT`                             | LAT drops from USFS. 
`MTBS`                            | Monitoring Trends in Burn Severity (MTBS) Burn Area Perimeters downloaded from: [USDA Forest Service FSGeodata Clearinghouse - Download National Datasets](https://data.fs.usda.gov/geodata/edw/datasets.php) under "MTBS Burn Area Boundary" (shapefile)
                                  | Ignition points (Fire Occurrence), and burn severity mosaics downloaded from [MTBS](https://www.mtbs.gov/direct-download).
`PARKS`                           | Any datasets taken from Parks, e.g. "DOBs_from_Parks"
`PARKS/DOBs_from_Parks`           | "DOBs_from_Parks" dataset
`QAQC`                            | 
`WFEIS`                           | 
`USFS`                            | USFS National forests `NationalForests`, roads `Roads`, wilderness areas `Wilderness` downloaded from [USDA Forest Service FSGeodata Clearinghouse - Download National Datasets](https://data.fs.usda.gov/geodata/edw/datasets.php) under "Administrative Forest Boundaries",  "National Forest System Roads", and "National Wilderness Areas".  
`Wen2023`                         | `clean` data folder taken from Wen et al. (2023) replication data. Downloaded from [Dropbox](https://www.dropbox.com/scl/fo/qb71jdvbc2r1zr8x22y6d/AKHgWH7t3MNfXMOc3L8pV3k?rlkey=iftnzzza6w1rqqw9yvnzdzg3w&e=1&dl=0)
`WUI`                             | Census block WUI polygons from Radeloff et al. (2023). Downloaded from [The 1990-2020 wildland-urban interface of the conterminous United States - geospatial data (4th Edition)](https://www.fs.usda.gov/rds/archive/catalog/RDS-2015-0012-4)
