# -*- coding: utf-8 -*-
"""
Script to plot Precipitation (Pr) from netCDF files for a specific date.
Clips seaward colors using inverted land polygon as a white overlay.
Set date in line 23

Created in April 2025
@author: SteynAS@ufs.ac.za
"""

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import matplotlib.patches as mpatches
from shapely.geometry import Polygon, MultiPolygon

# Set date string
YYYYMMDD = "20130119"

def load_pr_data(file_path, variable_name):
    # Load Pr data from NetCDF file
    dataset = nc.Dataset(file_path)
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    pr = dataset.variables[variable_name][:]
    return lon, lat, pr

def plot_pr(lon, lat, pr_data):
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})

    # Set the extent to match the data domain
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())

    # Add geographical features
    ax.add_feature(cfeature.LAND, edgecolor='black')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')

    # Define updated precipitation bins and colors
    bounds = [0, 0.1, 5, 10, 15, 25, 50, 100, 150, 200]
    colors = [
        '#ffffff',  # 0–0.1 mm (trace)
        '#e0f3db',  # 0.1–5 mm
        '#ccebc5',  # 5–10 mm
        '#a8ddb5',  # 10–15 mm
        '#41ab5d',  # dark green (15–25 mm)
        '#7bccc4',  # 25–50 mm
        '#4eb3d3',  # 50–100 mm
        '#2b8cbe',  # 100–150 mm
        '#0868ac',  # 150–200 mm
        '#191970'   # >200 mm - midnight blue (extend='max')
    ]

    cmap = mcolors.ListedColormap(colors)
    norm = mcolors.BoundaryNorm(bounds, len(bounds) - 1)

    # Mask only nodata values if needed
    pr_slice = pr_data[0, :, :]
    masked_pr = np.ma.masked_invalid(pr_slice)

    # Plot precipitation
    contour = ax.contourf(
        lon, lat, masked_pr,
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

    # Add colorbar with integer tick labels
    tick_labels = ['0.1', '5', '10', '15', '25', '50', '100', '150', '200']
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.05, ticks=bounds[1:])
    cbar.ax.set_yticklabels(tick_labels)
    cbar.set_label('Precipitation (mm/day)')

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
    pr_file = fr'C:\AgERA5\Pr\Precipitation-Flux_C3S-glob-agric_AgERA5_{YYYYMMDD}_final-v1.1.area-subset.-20.15.-35.35.nc'

    pr_variable = 'Precipitation_Flux'  # Updated variable name

    lon, lat, pr_data = load_pr_data(pr_file, pr_variable)

    # Plot Pr for specified date
    plot_pr(lon, lat, pr_data)
