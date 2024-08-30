# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 01:52:19 2024

@author: SteynAS
"""

import netCDF4 as nc
import matplotlib.pyplot as plt
import os
import glob
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def load_eto_data(file_path):
    dataset = nc.Dataset(file_path)
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    eto = dataset.variables['ETo'][:]  # 'ETo' is the variable name in the NetCDF file
    
    # Mask missing data
    eto = np.ma.masked_equal(eto, -9999)
    
    return lon, lat, eto

def process_files(directory, file_pattern):
    files = glob.glob(os.path.join(directory, file_pattern))
    
    if not files:
        raise FileNotFoundError(f"No files matching the pattern were found in the directory: {directory}")
    
    all_eto = []
    
    for file_path in files:
        lon, lat, eto = load_eto_data(file_path)
        all_eto.append(eto)
    
    # Concatenate all ETo data along a new axis (to simulate a time dimension)
    all_eto = np.ma.stack(all_eto, axis=0)
    
    return lon, lat, all_eto

def calculate_monthly_total(eto):
    # Sum daily ETo to get the total for the month (assuming each file is one day)
    monthly_total = np.ma.sum(eto, axis=0)
    return monthly_total

def plot_eto(lon, lat, eto, index=0):
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Set the extent to match the data domain
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())
    
    # Add geographical features
    ax.add_feature(cfeature.LAND, edgecolor='black')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')
    
    # Custom colour map
    cmap = plt.cm.viridis
    
    # Plot data
    eto_data = eto[index, :, :]  # Indexing over the new "time" dimension
    contour = ax.contourf(lon, lat, eto_data, cmap=cmap, transform=ccrs.PlateCarree())
    
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.05)
    cbar.set_label('ETo (mm/day)')
    
    ax.set_title(f'Reference Evapotranspiration (Day {index + 1})')
    
    # Add gridlines and labels
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    plt.show()

def plot_total_eto(lon, lat, total_eto):
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Set the extent to match the data domain
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())
    
    # Add geographical features
    ax.add_feature(cfeature.LAND, edgecolor='black')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')
    
    # Custom colour map
    cmap = plt.cm.viridis
    
    # Plot data
    contour = ax.contourf(lon, lat, total_eto, cmap=cmap, transform=ccrs.PlateCarree())
    
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.05)
    cbar.set_label('Total ETo (mm)')
    
    ax.set_title('Total Reference Evapotranspiration for January')
    
    # Add gridlines and labels
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    plt.show()

if __name__ == "__main__":
    directory = r'C:\AgERA5\ETo'  # Change to your ETo directory
    file_pattern = r'ETo_C3S-glob-agric_AgERA5_197901*_final-v1.1.area-subset.-20.15.-35.35.nc'
    
    lon, lat, all_eto = process_files(directory, file_pattern)
    
    # Calculate total ETo for January
    total_eto_january = calculate_monthly_total(all_eto)
    
    # Plot ETo for the first file (first day)
    plot_eto(lon, lat, all_eto, index=0)
    
    # Plot total ETo for January
    plot_total_eto(lon, lat, total_eto_january)
