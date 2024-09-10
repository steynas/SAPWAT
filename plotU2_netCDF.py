# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 03:05:51 2024

@author: SteynAS
"""

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def load_wind_speed_data(file_path, variable_name):
    dataset = nc.Dataset(file_path)
    time = dataset.variables['time'][:]
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    wind_speed = dataset.variables[variable_name][:]
    
    # Convert time to a readable format
    time_units = dataset.variables['time'].units
    time_calendar = dataset.variables['time'].calendar
    dates = nc.num2date(time, units=time_units, calendar=time_calendar)
    
    return lon, lat, dates, wind_speed

def reduce_wind_speed_to_2m_FAO56(wind_speed_10m):
    # FAO56 adjustment factor for converting 10m wind speed to 2m wind speed
    adjustment_factor = 4.87 / np.log(67.8 * 10 - 5.42)
    return wind_speed_10m * adjustment_factor

def plot_wind_speed(lon, lat, wind_speed, dates, index=0):
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Set the extent to match the data domain
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())
    
    # Add geographical features
    ax.add_feature(cfeature.LAND, edgecolor='black')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')
    
    # Custom colour map
    cmap = plt.get_cmap('viridis')
    
    # Plot data
    wind_speed_data = wind_speed[index, :, :]  # Wind speed is already in m/s
    contour = ax.contourf(lon, lat, wind_speed_data, cmap=cmap, transform=ccrs.PlateCarree())
    
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.05)
    cbar.set_label('Wind Speed at 2m (m/s)')
    
    ax.set_title(f'Wind Speed at 2m on {dates[index].strftime("%Y-%m-%d")}')
    
    # Add gridlines and labels
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    plt.show()

if __name__ == "__main__":
    # Path to the Wind Speed file
    ws_file = r'C:\AgERA5\WS\Wind-Speed-10m-Mean_C3S-glob-agric_AgERA5_19790101_final-v1.1.area-subset.-20.15.-35.35.nc'
    
    ws_variable = 'Wind_Speed_10m_Mean'
    
    lon, lat, dates_ws, wind_speed_10m = load_wind_speed_data(ws_file, ws_variable)
    
    # Reduce wind speed to 2m using the FAO56 formula
    wind_speed_2m = reduce_wind_speed_to_2m_FAO56(wind_speed_10m)
    
    # The expected date format in NetCDF after conversion should match the datetime format
    target_date = np.datetime64('1979-01-01')
    
    # Find the index for January 1st, 1979
    matching_indices = np.where(np.array([np.datetime64(date) for date in dates_ws]) == target_date)[0]
    
    if matching_indices.size == 0:
        print(f"Date {target_date} not found in the dataset.")
    else:
        target_index = matching_indices[0]
        # Plot WS at 2m for January 1st, 1979
        plot_wind_speed(lon, lat, wind_speed_2m, dates_ws, index=target_index)
