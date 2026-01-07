#### 01_download_wfeis_emissions.py ####
## Purpose: Downloads and extracts fire-level emissions data from WFEIS Calculator: https://wfeis.mtri.org/calculator

import os
import urllib.request
import pandas as pd
import geopandas as gpd
from io import StringIO

# Get the current directory (where the script is running)
current_dir = os.path.dirname(os.path.abspath(__file__))

# Navigate to the parent directory (LandscapeEffect)
project_dir = os.path.abspath(os.path.join(current_dir, '..', '..'))

# Now construct the path to the shapefile using the project root directory
mtbs_shapefile = os.path.join(project_dir, 'data', 'raw', 'MTBS', 'MTBS_Burn_Area', 'S_USA.MTBS_BURN_AREA_BOUNDARY.shp')

# Load the shapefile into a GeoDataFrame
gdf = gpd.read_file(mtbs_shapefile)

# Filter the data to include only burn areas from 2006 to 2023
filtered_gdf = gdf[(gdf['YEAR'] >= 2006) & (gdf['YEAR'] <= 2023)]

# Select only the FIRE_ID column
burn_ids = filtered_gdf['FIRE_ID'].tolist()

# Burned area data source
burnedarea_source = 'MTBS'

# Fuelbed system
fuelbed_system = 'fccs'

# Output format
output_format = 'csv'

# Output directory
# Now construct the path to the shapefile using the project root directory
output_dir = os.path.join(project_dir, 'data', 'raw', 'WFEIS')

# Create output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# List to hold all DataFrames
all_dfs = []

# Loop through each burn ID and download data
for burn_id in burn_ids:
    print(f'\nProcessing burn ID: {burn_id}')

    # Construct the URL for the WFEIS API request
    url = (
        f'https://wfeis.mtri.org/run_consume?'
        f'burn_ids={burn_id}&'
        f'burnedarea_source={burnedarea_source}&'
        f'shrub_blackened_pct=50&'
        f'fuelbed_system={fuelbed_system}&'
        f'outputs=fuelbed_aggregate&'
        f'output_format={output_format}&'
        f'output_spatial_aggregations=burnedarea&'
        f'output_temporal_aggregations=total'
    )

    print(f'URL: {url}')

    try:
        # Send the request and read the response
        with urllib.request.urlopen(url) as response:
            data = response.read().decode('utf-8')

        # Convert the response data to a pandas DataFrame
        df = pd.read_csv(StringIO(data))

        # Add a new column 'FIRE_ID' with the current burn ID
        df['FIRE_ID'] = burn_id

        # Rename the columns to match the required names
        df.rename(columns={'consume_output__pm25_mg': 'FIRE_PM25', 'consume_output__co2_mg': 'FIRE_CO2'}, inplace=True)

        # Select only the required columns: FIRE_ID, FIRE_PM25, and FIRE_CO2
        df = df[['FIRE_ID', 'FIRE_PM25', 'FIRE_CO2']]

        # Append the DataFrame to the list
        all_dfs.append(df)

        print(f'Data successfully processed for burn ID: {burn_id}')
    except Exception as e:
        print(f'Error retrieving data for burn ID {burn_id}: {e}')

# Combine all DataFrames into one
combined_df = pd.concat(all_dfs, ignore_index=True)

# Save the combined DataFrame to a single CSV file
combined_output_file = os.path.join(output_dir, 'WFEIS_data.csv')
combined_df.to_csv(combined_output_file, index=False)

print(f'All data successfully combined and written to {combined_output_file}')

