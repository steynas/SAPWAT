# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 02:57:28 2024

@author: SteynAS
"""

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def load_temperature_data(file_path, variable_name):
    dataset = nc.Dataset(file_path)
    time = dataset.variables['time'][:]
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    temperature = dataset.variables[variable_name][:]
    
    # Convert time to a readable format
    time_units = dataset.variables['time'].units
    time_calendar = dataset.variables['time'].calendar
    dates = nc.num2date(time, units=time_units, calendar=time_calendar)
    
    return lon, lat, dates, temperature

def plot_temperature(lon, lat, temperature, dates, index=0):
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Set the extent to match the data domain
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())
    
    # Add geographical features
    ax.add_feature(cfeature.LAND, edgecolor='black')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')
    
    # Custom colour map
    cmap = plt.get_cmap('coolwarm')
    
    # Plot data
    temperature_data = temperature[index, :, :] - 273.15  # Convert from Kelvin to Celsius
    contour = ax.contourf(lon, lat, temperature_data, cmap=cmap, transform=ccrs.PlateCarree())
    
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.05)
    cbar.set_label('Mean Temperature (°C)')
    
    ax.set_title(f'Mean Temperature on {dates[index].strftime("%Y-%m-%d")}')
    
    # Add gridlines and labels
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    plt.show()

if __name__ == "__main__":
    # Updated paths with the correct directory structure
    tmin_file = r'C:\AgERA5\Temps\Temperature-Air-2m-Min-24h_C3S-glob-agric_AgERA5_19790101_final-v1.1.area-subset.-20.15.-35.35.nc'
    tmax_file = r'C:\AgERA5\Temps\Temperature-Air-2m-Max-24h_C3S-glob-agric_AgERA5_19790101_final-v1.1.area-subset.-20.15.-35.35.nc'
    
    tmin_variable = 'Temperature_Air_2m_Min_24h'
    tmax_variable = 'Temperature_Air_2m_Max_24h'
    
    lon, lat, dates_tmin, tmin = load_temperature_data(tmin_file, tmin_variable)
    _, _, dates_tmax, tmax = load_temperature_data(tmax_file, tmax_variable)
    
    # Assuming the dates in Tmin and Tmax files match, calculate Tmean
    tmean = (tmax + tmin) / 2
    
    # The expected date format in NetCDF after conversion should match the datetime format
    target_date = np.datetime64('1979-01-01')
    
    # Find the index for January 1st, 1979
    matching_indices = np.where(np.array([np.datetime64(date) for date in dates_tmin]) == target_date)[0]
    
    if matching_indices.size == 0:
        print(f"Date {target_date} not found in the dataset.")
    else:
        target_index = matching_indices[0]
        # Plot Tmean for January 1st, 1979
        plot_temperature(lon, lat, tmean, dates_tmin, index=target_index)
