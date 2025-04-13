# -*- coding: utf-8 -*-
"""
Script to plot Mean Vapour Pressure (VP) from netCDF files for a specific date.
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

def load_vp_data(file_path, variable_name):
    dataset = nc.Dataset(file_path)
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    vp = dataset.variables[variable_name][:]  # (time, lat, lon)
    return lon, lat, vp

def plot_vp(lon, lat, vp_data):
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.LAND, edgecolor='black')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')

    # Define bins and colours (VP >= 1.0 hPa), final bin ends at 36, extend='max'
    bounds = [0, 1, 4, 8, 12, 16, 20, 24, 28, 32, 36]
    colours = [
        '#ffffff',  # VP < 1.0 (white)
        '#deebf7', '#c6dbef', '#9ecae1', '#6baed6',
        '#4292c6', '#2171b5', '#08519c', '#084594',  # 28–32
        '#08306b',  # 32–36
        '#000000'   # >36 (extend='max')
    ]

    cmap = mcolors.ListedColormap(colours)
    norm = mcolors.BoundaryNorm(bounds, len(bounds) - 1)

    data = vp_data[0, :, :]
    masked_vp = np.ma.masked_invalid(data)

    # Plot all bins including VP < 1.0
    contour = ax.contourf(
        lon, lat, masked_vp,
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

    tick_labels = ['1', '4', '8', '12', '16', '20', '24', '28', '32', '36']
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.05,
                        ticks=bounds[1:])
    cbar.ax.set_yticklabels(tick_labels)
    cbar.set_label('Mean Vapour Pressure (hPa)')

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlocator = mticker.MultipleLocator(2)
    gl.ylocator = mticker.MultipleLocator(2)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}

    plt.show()

if __name__ == "__main__":
    vp_file = fr'C:\AgERA5\VP\Vapour-Pressure-Mean_C3S-glob-agric_AgERA5_{YYYYMMDD}_final-v1.1.area-subset.-20.15.-35.35.nc'
    vp_variable = 'Vapour_Pressure_Mean'

    lon, lat, vp_data = load_vp_data(vp_file, vp_variable)

    # Plot Mean Vapour Pressure for specified date
    plot_vp(lon, lat, vp_data)