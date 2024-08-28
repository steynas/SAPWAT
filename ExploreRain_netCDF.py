# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 13:34:10 2024

@author: SteynAS
"""

import xarray as xr

# Path to the Precipitation Flux NetCDF file
file_path = r"C:\AgERA5\Rain\Precipitation-Flux_C3S-glob-agric_AgERA5_19790101_final-v1.1.area-subset.-20.15.-35.35.nc"

# Open the NetCDF file
try:
    ds = xr.open_dataset(file_path)
    print("Dataset successfully opened.")
    print(ds)
except Exception as e:
    print(f"An error occurred: {e}")

# Print out the dataset details to see variable names and structure
print(ds)
print(ds['time'].values)