# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 14:52:01 2024

@author: SteynAS
"""

import rasterio
import pandas as pd
import os

# Define file paths
base_dir = "C:/AgERA5"
dem_dir = os.path.join(base_dir, "DEM")
dem_file = os.path.join(dem_dir, "merged_dem.tif")
grid_points_file = os.path.join(dem_dir, "AgERA5_GridPoints.csv")
output_file = os.path.join(dem_dir, "AgERA5_GridPoint_Elevations.csv")

# Load the DEM file
with rasterio.open(dem_file) as dem:
    dem_data = dem.read(1)  # Load the first band (elevation data)
    dem_transform = dem.transform
    dem_bounds = dem.bounds  # Get the bounds of the DEM

# Load the grid points CSV file
grid_points = pd.read_csv(grid_points_file)

# Function to get elevation at a given lat/lon and force replace -32768 with "NaN"
def get_elevation(lat, lon, dem_data, dem_transform, dem_bounds):
    # Check if the point is within the bounds of the DEM
    if not (dem_bounds.left <= lon <= dem_bounds.right and dem_bounds.bottom <= lat <= dem_bounds.top):
        print(f"Point ({lat}, {lon}) is out of DEM bounds.")
        return "NaN"
    
    lat = round(lat, 4)  # Use higher precision
    lon = round(lon, 4)
    
    try:
        row, col = rasterio.transform.rowcol(dem_transform, lon, lat)
        elevation = dem_data[row, col]
        if elevation == -32768:
            return "NaN"
        return elevation
    except IndexError:
        # This happens if the coordinates are slightly out of the DEM grid
        print(f"Point ({lat}, {lon}) is out of DEM grid.")
        return "NaN"

# Apply the function to all grid points
grid_points['elevation'] = [
    get_elevation(lat, lon, dem_data, dem_transform, dem_bounds) for lat, lon in zip(grid_points['latitude'], grid_points['longitude'])
]

# Save the results to a new CSV file
grid_points.to_csv(output_file, index=False, na_rep="NaN")

print(f'Elevation data saved to {output_file}')

