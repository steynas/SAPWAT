# -*- coding: utf-8 -*-
"""
Script to plot Reference Evapotranspiration (ETo) for 1979-01-01
"""

import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def load_eto_data(file_path, variable_name):
    # Load ETo data from NetCDF file
    dataset = nc.Dataset(file_path)
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    eto = dataset.variables[variable_name][:]

    return lon, lat, eto

def plot_eto(lon, lat, eto_data):
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Set the extent to match the data domain
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())
    
    # Add geographical features
    ax.add_feature(cfeature.LAND, edgecolor='black')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')
    
    # Custom color map
    cmap = plt.get_cmap('viridis')
    
    # Plot ETo data (assuming it's 3D with time dimension first)
    contour = ax.contourf(lon, lat, eto_data[0, :, :], cmap=cmap, transform=ccrs.PlateCarree())
    
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.05)
    cbar.set_label('Reference Evapotranspiration (mm/day)')
    
    # Add gridlines and labels
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    
    plt.show()

if __name__ == "__main__":
    # Path to the ETo file
    eto_file = r'C:\AgERA5\ETo\ETo_C3S-glob-agric_AgERA5_19790101_final-v1.1.area-subset.-20.15.-35.35.nc'
    
    eto_variable = 'Reference_Evapotranspiration_(PM-FAO56)'  # Updated variable name
    
    lon, lat, eto_data = load_eto_data(eto_file, eto_variable)
    
    # Plot ETo for 1979-01-01 (assuming it's 3D with time as the first dimension)
    plot_eto(lon, lat, eto_data)
