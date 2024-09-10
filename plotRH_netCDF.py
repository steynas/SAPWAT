# -*- coding: utf-8 -*-
"""
Adapted to plot RHmin and RHmax for 1979-01-01
"""

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Function to load RH data from the netCDF file
def load_rh_data(file_path, variable_name):
    dataset = nc.Dataset(file_path)
    time = dataset.variables['time'][:]
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    rh_data = dataset.variables[variable_name][:]  # Load the RH data
    
    # Convert time to a readable format
    time_units = dataset.variables['time'].units
    time_calendar = dataset.variables['time'].calendar
    dates = nc.num2date(time, units=time_units, calendar=time_calendar)
    
    return lon, lat, dates, rh_data

# Function to plot the RH data
def plot_rh(lon, lat, rh_data, dates, index=0, title="Relative Humidity"):
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
    
    # Plot RH data
    rh_data_on_date = rh_data[index, :, :]  # RH on the target date
    contour = ax.contourf(lon, lat, rh_data_on_date, cmap=cmap, transform=ccrs.PlateCarree())
    
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.05)
    cbar.set_label('Relative Humidity (%)')
    
    ax.set_title(f'{title} on {dates[index].strftime("%Y-%m-%d")}')
    
    # Add gridlines and labels
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    plt.show()

if __name__ == "__main__":
    # Path to the RHmin and RHmax files
    directory = r"C:/AgERA5/RH"
    rhmin_file = f"{directory}/RHmin-2m_C3S-glob-agric_AgERA5_19790101_final-v1.1.area-subset.-20.15.-35.35.nc"
    rhmax_file = f"{directory}/RHmax-2m_C3S-glob-agric_AgERA5_19790101_final-v1.1.area-subset.-20.15.-35.35.nc"
    
    # Load RHmin
    rhmin_variable = 'RHmin'  # Assuming 'RHmin' is the variable name in the file
    lon, lat, dates_rhmin, rhmin_data = load_rh_data(rhmin_file, rhmin_variable)
    
    # Load RHmax
    rhmax_variable = 'RHmax'  # Assuming 'RHmax' is the variable name in the file
    _, _, dates_rhmax, rhmax_data = load_rh_data(rhmax_file, rhmax_variable)
    
    # The expected date format in NetCDF after conversion should match the datetime format
    target_date = np.datetime64('1979-01-01')
    
    # Find the index for January 1st, 1979 in RHmin and RHmax files
    matching_indices_rhmin = np.where(np.array([np.datetime64(date) for date in dates_rhmin]) == target_date)[0]
    matching_indices_rhmax = np.where(np.array([np.datetime64(date) for date in dates_rhmax]) == target_date)[0]
    
    if matching_indices_rhmin.size == 0:
        print(f"Date {target_date} not found in the RHmin dataset.")
    else:
        target_index_rhmin = matching_indices_rhmin[0]
        # Plot RHmin for January 1st, 1979
        plot_rh(lon, lat, rhmin_data, dates_rhmin, index=target_index_rhmin, title="RHmin")
    
    if matching_indices_rhmax.size == 0:
        print(f"Date {target_date} not found in the RHmax dataset.")
    else:
        target_index_rhmax = matching_indices_rhmax[0]
        # Plot RHmax for January 1st, 1979
        plot_rh(lon, lat, rhmax_data, dates_rhmax, index=target_index_rhmax, title="RHmax")
