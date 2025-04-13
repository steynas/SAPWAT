# -*- coding: utf-8 -*-
"""
Script to calculate and plot Mean Temperature (Tmean) from netCDF files 
containing Tmin and Tmax for a specific date.
Clips seaward colors using inverted land polygon as a white overlay.
Set date in line 24

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

# Define bin edges for 9 bins between -10 and 35
bounds = [-10, -5, 0, 5, 10, 15, 20, 25, 30, 35]

# Define 9 bin colors only (excluding under/over colors)
bin_colors = [
    '#2171b5', '#4292c6', '#6baed6',
    '#9ecae1', '#c6dbef', '#fdd0a2',
    '#fdae6b', '#f16913', '#d94801'
]

# Create colormap and explicitly assign under/over colors
cmap = mcolors.ListedColormap(bin_colors)
cmap.set_under('#08306b')   # < -10°C
cmap.set_over('#67000d')    # > 35°C
norm = mcolors.BoundaryNorm(bounds, ncolors=len(bin_colors))

def load_temperature_data(file_path, variable_name):
    dataset = nc.Dataset(file_path)
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    temp = dataset.variables[variable_name][:]

    temp_2d = temp[0, :, :]  # First timestep
    temp_2d_celsius = temp_2d - 273.15
    return lon, lat, temp_2d_celsius

def plot_temperature(lon, lat, temp_data):
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())

    # Add map features
    ax.add_feature(cfeature.LAND, facecolor='none', edgecolor='black', zorder=4)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linestyle=':', zorder=5)
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none', zorder=5)
    ax.add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=0.8, zorder=6)

    # Force triangle overflow to appear for validation
    temp_data[0, 0] = -11.0
    temp_data[0, 1] = 36.0

    # Masked temperature data
    masked_temp = np.ma.masked_invalid(temp_data)

    # Plot with contourf
    contour = ax.contourf(
        lon, lat, masked_temp,
        levels=bounds, cmap=cmap, norm=norm,
        extend='both', transform=ccrs.PlateCarree(),
        zorder=3
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
                facecolor='white', edgecolor='none', transform=ccrs.PlateCarree(), zorder=4
            )
            ax.add_patch(patch)

    # Add colorbar
    tick_values = bounds
    tick_labels = [str(t) for t in tick_values]
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.05, extend='both', ticks=tick_values)
    cbar.set_label('Temperature (°C)')
    cbar.ax.set_yticklabels(tick_labels)

    # Gridlines
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlocator = mticker.MultipleLocator(2)
    gl.ylocator = mticker.MultipleLocator(2)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}

    # Add black box around the entire plot area
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.2)

    plt.show()

if __name__ == "__main__":
    tmin_file = fr'C:\\AgERA5\\Temps\\Temperature-Air-2m-Min-24h_C3S-glob-agric_AgERA5_{YYYYMMDD}_final-v1.1.area-subset.-20.15.-35.35.nc'
    tmax_file = fr'C:\\AgERA5\\Temps\\Temperature-Air-2m-Max-24h_C3S-glob-agric_AgERA5_{YYYYMMDD}_final-v1.1.area-subset.-20.15.-35.35.nc'

    tmin_variable = 'Temperature_Air_2m_Min_24h'
    tmax_variable = 'Temperature_Air_2m_Max_24h'

    lon_tmin, lat_tmin, tmin_data = load_temperature_data(tmin_file, tmin_variable)
    lon_tmax, lat_tmax, tmax_data = load_temperature_data(tmax_file, tmax_variable)

    tmean_data = (tmin_data + tmax_data) / 2.0
    plot_temperature(lon_tmin, lat_tmin, tmean_data)

