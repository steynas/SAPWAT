# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 06:23:05 2024

@author: SteynAS
"""

import xarray as xr

# Path to the NetCDF file
file_path = 'C:\AgERA5\Temperature\Temperature-Air-2m-Max-24h_C3S-glob-agric_AgERA5_19790101_final-v1.1.area-subset.-20.15.-35.35.nc'

# Open the NetCDF file
ds = xr.open_dataset(file_path)

# Print out the dataset details to see variable names and structure
print(ds)
