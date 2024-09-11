# -*- coding: utf-8 -*-
"""
This is a script to extract AgERA5 grid point coordinates
and writhe them to a CSV file

@author: SteynAS@ufs.ac.za
"""

import netCDF4 as nc
import csv
import os

# Set base directory
base_dir = "C:/AgERA5"
vp_dir = os.path.join(base_dir, "VP")
dem_dir = os.path.join(base_dir, "DEM")
output_csv_path = os.path.join(dem_dir, "AgERA5_GridPoints.csv")

# Path to the NetCDF file
nc_file_path = os.path.join(vp_dir, "Vapour-Pressure-Mean_C3S-glob-agric_AgERA5_19790101_final-v1.1.area-subset.-20.15.-35.35.nc")

# Open the NetCDF file
dataset = nc.Dataset(nc_file_path)

# Extract latitude and longitude data from the NetCDF file
latitudes = dataset.variables['lat'][:]
longitudes = dataset.variables['lon'][:]

# Writing latitudes and longitudes to the CSV file
with open(output_csv_path, mode='w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(["latitude", "longitude"])  # Writing the header

    for lat in latitudes:
        for lon in longitudes:
            writer.writerow([lat, lon])

print(f"Grid points have been written to {output_csv_path}")

# Close the NetCDF file
dataset.close()
