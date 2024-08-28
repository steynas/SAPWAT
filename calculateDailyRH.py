# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 17:57:17 2024

@author: SteynAS
"""

import xarray as xr
import numpy as np
import os
import calendar

# Request user input for year and month
year = input("Enter the year (e.g., 1979): ")
month = input("Enter the month (01-12): ")

# Determine if the selected year is a leap year
is_leap_year = calendar.isleap(int(year))

# Determine the number of days in the selected month
days_in_month = calendar.monthrange(int(year), int(month))[1]

# Directory containing the NetCDF files
directory = 'C:/AgERA5/RH/'

# Define the study area coordinates
min_lat, max_lat = -35, -20
min_lon, max_lon = 15, 35

# List of potential variable names
time_to_variable = {
    '06:00': 'Relative_Humidity_2m_06h',
    '09:00': 'Relative_Humidity_2m_09h',
    '12:00': 'Relative_Humidity_2m_12h',
    '15:00': 'Relative_Humidity_2m_15h',
    '18:00': 'Relative_Humidity_2m_18h'
}

# Iterate over each day in the selected month
for day in range(1, days_in_month + 1):
    target_date = f"{year}{month.zfill(2)}{str(day).zfill(2)}"

    # Initialize dictionaries to hold the RH values for all grid points
    rh_min_values = None
    rh_max_values = None

    # Iterate over each file in the directory
    for filename in sorted(os.listdir(directory)):
        if filename.endswith(".nc") and f"_{target_date}_" in filename:
            file_path = os.path.join(directory, filename)
            
            # Open the NetCDF file
            ds = xr.open_dataset(file_path)
            
            # Check which time this file corresponds to
            for time, variable_name in time_to_variable.items():
                if variable_name in ds:
                    print(f"Processing {variable_name} in {filename}")
                    
                    # Use the correct names for latitude and longitude
                    lats = ds['lat'].values
                    lons = ds['lon'].values
                    
                    # Determine the indices within the study area
                    lat_indices = np.where((lats >= min_lat) & (lats <= max_lat))[0]
                    lon_indices = np.where((lons >= min_lon) & (lons <= max_lon))[0]
                    
                    # Extract RH data for the study area
                    rh_data = ds[variable_name].values[0, lat_indices.min():lat_indices.max()+1, lon_indices.min():lon_indices.max()+1]
                    
                    # Initialize min and max arrays on the first run
                    if rh_min_values is None:
                        rh_min_values = np.full(rh_data.shape, np.inf)
                        rh_max_values = np.full(rh_data.shape, -np.inf)
                    
                    # Update min and max values, ignoring NaNs
                    rh_min_values = np.fmin(rh_min_values, np.nan_to_num(rh_data, nan=np.inf))
                    rh_max_values = np.fmax(rh_max_values, np.nan_to_num(rh_data, nan=-np.inf))
                    
                    break

    # Replace inf values with NaN to reflect areas with no valid data
    rh_min_values[rh_min_values == np.inf] = np.nan
    rh_max_values[rh_max_values == -np.inf] = np.nan

    # Create a template for the new .nc files using one of the datasets
    output_ds_min = xr.Dataset(
        {
            "RHmin": (["lat", "lon"], rh_min_values),
        },
        coords={
            "lat": ds['lat'].values[lat_indices.min():lat_indices.max()+1],
            "lon": ds['lon'].values[lon_indices.min():lon_indices.max()+1]
        },
        attrs=ds.attrs  # Copy global attributes from the input dataset
    )

    output_ds_max = xr.Dataset(
        {
            "RHmax": (["lat", "lon"], rh_max_values),
        },
        coords={
            "lat": ds['lat'].values[lat_indices.min():lat_indices.max()+1],
            "lon": ds['lon'].values[lon_indices.min():lon_indices.max()+1]
        },
        attrs=ds.attrs  # Copy global attributes from the input dataset
    )

    # Define the output directory (same as input directory)
    output_dir = directory

    # Generate filenames based on the input filename pattern and save them for the specified date
    rhmin_filename = os.path.join(output_dir, f"RHmin-2m_C3S-glob-agric_AgERA5_{target_date}_final-v1.1.area-subset.-20.15.-35.35.nc")
    rhmax_filename = os.path.join(output_dir, f"RHmax-2m_C3S-glob-agric_AgERA5_{target_date}_final-v1.1.area-subset.-20.15.-35.35.nc")

    # Write out the new .nc files using the "netcdf4" engine to avoid the cfgrib issue
    output_ds_min.to_netcdf(rhmin_filename, engine="netcdf4")
    output_ds_max.to_netcdf(rhmax_filename, engine="netcdf4")

    print(f"RHmin saved to {rhmin_filename}")
    print(f"RHmax saved to {rhmax_filename}")
