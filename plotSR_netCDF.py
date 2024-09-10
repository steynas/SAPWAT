# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 03:12:52 2024

@author: SteynAS
"""

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def load_solar_radiation_data(file_path, variable_name):
    dataset = nc.Dataset(file_path)
    time = dataset.variables['time'][:]
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    solar_radiation_flux = dataset.variables[variable_name][:]
    
    # Convert time to a readable format
    time_units = dataset.variables['time'].units
    time_calendar = dataset.variables['time'].calendar
    dates = nc.num2date(time, units=time_units, calendar=time_calendar)
    
    return lon, lat, dates, solar_radiation_flux

def plot_solar_radiation(lon, lat, solar_radiation, dates, index=0):
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Set the extent to match the data domain
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())
    
    # Add geographical features
    ax.add_feature(cfeature.LAND, edgecolor='black')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')
    
    # Custom colour map
    cmap = plt.get_cmap('plasma')
    
    # Plot data
    solar_radiation_data = solar_radiation[index, :, :]  # Original units (J/m²/day)
    contour = ax.contourf(lon, lat, solar_radiation_data, cmap=cmap, transform=ccrs.PlateCarree())
    
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.05)
    cbar.set_label('Solar Radiation (J/m²/day)')
    
    # Remove the title (header)
    # ax.set_title(f'Solar Radiation on {dates[index].strftime("%Y-%m-%d")}')
    
    # Add gridlines and labels
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    plt.show()

if __name__ == "__main__":
    # Path to the Solar Radiation file
    sr_file = r'C:\AgERA5\SR\Solar-Radiation-Flux_C3S-glob-agric_AgERA5_19790101_final-v1.1.area-subset.-20.15.-35.35.nc'
    
    sr_variable = 'Solar_Radiation_Flux'
    
    lon, lat, dates_sr, solar_radiation_flux = load_solar_radiation_data(sr_file, sr_variable)
    
    # The expected date format in NetCDF after conversion should match the datetime format
    target_date = np.datetime64('1979-01-01')
    
    # Find the index for January 1st, 1979
    matching_indices = np.where(np.array([np.datetime64(date) for date in dates_sr]) == target_date)[0]
    
    if matching_indices.size == 0:
        print(f"Date {target_date} not found in the dataset.")
    else:
        target_index = matching_indices[0]
        # Plot SR for January 1st, 1979 in original units (J/m²/day)
        plot_solar_radiation(lon, lat, solar_radiation_flux, dates_sr, index=target_index)
