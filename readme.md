# Data Sources

## Setting Up "data/raw" directory

In order for the code to run we need to have the same folder/file names in "data/raw" for the data downloaded from the same links: 

Data folder     | Description
--------------------|------------------------------------------
`ACS`                             | American Community Survey (ACS) Census Block Group Housing Data downloaded from: [IPUMS](https://usa.ipums.org/usa/)
`FACTS`                           | USFS Fuel Treatment Data downloaded from: [USDA Forest Service FSGeodata Clearinghouse - Download National Datasets]https://data.fs.usda.gov/geodata/edw/datasets.php under "Hazardous Fuel Treatment Reduction: Polygon" (shapefile) - Rename downloaded folder to "FACTS"
`MTBS`                            | MTBS Burn Area Perimeters downloaded from: [USDA Forest Service FSGeodata Clearinghouse - Download National Datasets]https://data.fs.usda.gov/geodata/edw/datasets.php under "MTBS Burn Area Boundary" (shapefile) - Rename downloaded folder to "MTBS"
`WUI`                             | Census block WUI polygons downloaded from: https://www.fs.usda.gov/rds/archive/catalog/RDS-2015-0012-4 - Rename downloaded folder to "WUI"
`PARKS`                           | Any datasets taken from Parks, e.g. "DOBs_from_Parks"
`PARKS/DOBs_from_Parks`           | "DOBs_from_Parks" dataset
`hourly_perims`                   | hourly wildfire perimeters dataset
`FIRED`                           | FIRED daily wildfire perimeters dataset downloaded from: https://scholar.colorado.edu/concern/datasets/8336h304x
`Highways`                        | U.S. Highway roads downloaded from: https://catalog.data.gov/dataset/tiger-line-shapefile-2016-nation-u-s-primary-roads-national-shapefile under Shapefile Zip File
`USFS_Roads`                      | USFS Roads downloaded from: https://data.fs.usda.gov/geodata/edw/datasets.php under "National Forest System Roads" 
`USFS_Boundaries/NationalForests` | USFS National Forest Boundaries downloaded from: https://data.fs.usda.gov/geodata/edw/datasets.php under "Administrative Forest Boundaries" 
`USFS_Boundaries/WildernessAreas` | USFS Wilderness Areas downloaded from: https://data.fs.usda.gov/geodata/edw/datasets.php under "National Wilderness Areas" 
