# -*- coding: utf-8 -*-
"""
Script to plot Minimum Temperature (Tmin) from netCDF files for a specific date.
Clips seaward colors using inverted land polygon as a white overlay.
Rogue Tmin values <-30°C can be validated visually.
Set date in line 25

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
from shapely.geometry import Polygon
from shapely.ops import unary_union

# Set target date string
YYYYMMDD = "20231127"

# Define bin edges and colours for Tmin
bounds = [-30, -20, -10, -5, 0, 5, 10, 15, 20, 25, 30, 35]
bin_colors = [
    '#08306b', '#2171b5', '#4292c6', '#6baed6',
    '#9ecae1', '#c6dbef', '#fdd0a2', '#fdae6b',
    '#f16913', '#d94801', '#a63603', '#67000d'
]
cmap = mcolors.ListedColormap(bin_colors)
cmap.set_under('#000000')  # below -30°C
cmap.set_over('#4d004b')   # above 35°C
norm = mcolors.BoundaryNorm(bounds, ncolors=len(bin_colors))

def load_temperature_data(file_path, variable_name):
    dataset = nc.Dataset(file_path)
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    temp = dataset.variables[variable_name][:]
    temp_2d = temp[0, :, :] - 273.15  # First timestep, convert to Celsius
    return lon, lat, temp_2d

def plot_temperature(lon, lat, temp_data):
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())

    # Create meshgrid extent polygon for ocean masking
    lon2d, lat2d = np.meshgrid(lon, lat)
    lon_min, lon_max = lon2d.min(), lon2d.max()
    lat_min, lat_max = lat2d.min(), lat2d.max()
    full_extent = Polygon([
        (lon_min, lat_min), (lon_max, lat_min),
        (lon_max, lat_max), (lon_min, lat_max),
        (lon_min, lat_min)
    ])

    # Land union and mask creation
    land_shp = shpreader.natural_earth(resolution='10m', category='physical', name='land')
    land_geoms = list(shpreader.Reader(land_shp).geometries())
    land_union = unary_union(land_geoms)
    ocean_geom = full_extent.difference(land_union)

    # Plot temperature
    masked_temp = np.ma.masked_invalid(temp_data)
    contour = ax.contourf(
        lon, lat, masked_temp,
        levels=bounds, cmap=cmap, norm=norm,
        extend='both', transform=ccrs.PlateCarree(), zorder=2
    )

    # Add white ocean patch
    if ocean_geom.geom_type == 'Polygon':
        patch = mpatches.Polygon(
            list(ocean_geom.exterior.coords), facecolor='white', edgecolor='none',
            transform=ccrs.PlateCarree(), zorder=3
        )
        ax.add_patch(patch)
    elif ocean_geom.geom_type == 'MultiPolygon':
        for poly in ocean_geom.geoms:
            patch = mpatches.Polygon(
                list(poly.exterior.coords), facecolor='white', edgecolor='none',
                transform=ccrs.PlateCarree(), zorder=3
            )
            ax.add_patch(patch)

    # Add map features on top
    ax.add_feature(cfeature.LAND, facecolor='none', edgecolor='black', zorder=4)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linestyle=':', zorder=5)
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none', zorder=5)
    ax.add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=0.8, zorder=6)

    # Add colourbar
    tick_labels = [str(t) for t in bounds]
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.05, extend='both', ticks=bounds)
    cbar.set_label('Minimum Temperature (°C)')
    cbar.ax.set_yticklabels(tick_labels)

    # Gridlines
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlocator = mticker.MultipleLocator(2)
    gl.ylocator = mticker.MultipleLocator(2)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}

    # Black box around plot
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.2)

    plt.title(f'AgERA5 Tmin on {YYYYMMDD}')
    plt.show()

if __name__ == "__main__":
    tmin_file = fr'C:\\AgERA5\\Temps\\Temperature-Air-2m-Min-24h_C3S-glob-agric_AgERA5_{YYYYMMDD}_final-v1.1.area-subset.-20.15.-35.35.nc'
    tmin_variable = 'Temperature_Air_2m_Min_24h'
    lon, lat, tmin_data = load_temperature_data(tmin_file, tmin_variable)
    plot_temperature(lon, lat, tmin_data)
