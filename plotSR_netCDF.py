# -*- coding: utf-8 -*-
"""
Script to plot Solar Radiation (SR) from netCDF files for a specific date.
Clips seaward colors using inverted land polygon as a white overlay.
Set date in line 23

Created in April 2025
@author: SteynAS@ufs.ac.za
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import cartopy.io.shapereader as shpreader
import matplotlib.patches as mpatches
from shapely.geometry import Polygon, MultiPolygon

# Set date string
YYYYMMDD = "20130119"

# Define SR bins and colours in MJ/m²/day (converted from J by /1e6)
bounds = [0, 5, 10, 15, 20, 25, 30, 35]

# Get a subset of the plasma colormap for bins and adjust colors for brighter last bin
plasma_cmap = plt.cm.plasma
colors = [plasma_cmap(i) for i in np.linspace(0, 0.85, len(bounds)-1)]  # Adjusted to not go to the darkest colour

# Use the brightest colour in plasma for >35
colors.append('#ffff00')  # Bright yellow colour for >35 MJ/m²/day (most distinct)

cmap = mcolors.ListedColormap(colors)
norm = mcolors.BoundaryNorm(bounds, len(bounds) - 1)

def load_sr_data(file_path, variable_name):
    dataset = nc.Dataset(file_path)
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    sr = dataset.variables[variable_name][:]
    return lon, lat, sr

def plot_sr(lon, lat, sr_data):
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})

    # Set the extent to match the data domain
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())

    # Add geographical features
    ax.add_feature(cfeature.LAND, edgecolor='black')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')

    # Convert SR from J/m²/day to MJ/m²/day
    sr_mj = sr_data[0, :, :] / 1e6  # Conversion to MJ/m²/day
    masked_sr = np.ma.masked_invalid(sr_mj)

    # Plot solar radiation data
    contour = ax.contourf(
        lon, lat, masked_sr,
        levels=bounds,
        cmap=cmap, norm=norm, extend='max', transform=ccrs.PlateCarree()
    )


    # White ocean mask (using inverted land polygons)
    land_shp = shpreader.natural_earth(resolution='10m', category='physical', name='land')
    geometries = list(shpreader.Reader(land_shp).geometries())

    # Define full plotting extent as outer box
    full_extent = Polygon([
        (lon.min(), lat.min()),
        (lon.max(), lat.min()),
        (lon.max(), lat.max()),
        (lon.min(), lat.max()),
        (lon.min(), lat.min())
    ])

    # Subtract all land geometries from full extent
    if len(geometries) == 1:
        land_union = geometries[0]
    else:
        land_parts = []
        for g in geometries:
            if g.geom_type == 'Polygon':
                land_parts.append(g)
            elif g.geom_type == 'MultiPolygon':
                land_parts.extend([poly for poly in g.geoms])
        land_union = MultiPolygon(land_parts)

    ocean_geom = full_extent.difference(land_union)

    # Add white ocean patch
    if ocean_geom.geom_type == 'Polygon':
        patch = mpatches.Polygon(
            list(ocean_geom.exterior.coords),
            facecolor='white', edgecolor='none', transform=ccrs.PlateCarree(), zorder=1
        )
        ax.add_patch(patch)
    elif ocean_geom.geom_type == 'MultiPolygon':
        for poly in ocean_geom.geoms:
            patch = mpatches.Polygon(
                list(poly.exterior.coords),
                facecolor='white', edgecolor='none', transform=ccrs.PlateCarree(), zorder=1
            )
            ax.add_patch(patch)

    # Add colorbar with appropriate ticks
    tick_labels = ['0', '5', '10', '15', '20', '25', '30', '35']
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.05, ticks=bounds)
    cbar.ax.set_yticklabels(tick_labels)
    cbar.set_label('Solar Radiation (MJ/m²/day)')

    # Add gridlines and labels with 2x2 degree spacing
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlocator = mticker.MultipleLocator(2)
    gl.ylocator = mticker.MultipleLocator(2)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}

    plt.show()

if __name__ == "__main__":
    # Construct file path from date
    sr_file = fr'C:\AgERA5\SR\Solar-Radiation-Flux_C3S-glob-agric_AgERA5_{YYYYMMDD}_final-v1.1.area-subset.-20.15.-35.35.nc'
    sr_variable = 'Solar_Radiation_Flux'  # Variable name for SR
    
    lon, lat, sr_data = load_sr_data(sr_file, sr_variable)

    # Plot SR for specified date
    plot_sr(lon, lat, sr_data)
