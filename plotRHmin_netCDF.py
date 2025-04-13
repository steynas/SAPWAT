# -*- coding: utf-8 -*-
"""
Script to plot Minimum Relative Humidity (RHmin) from netCDF files for a specific date.
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

def load_rhmin_data(file_path, variable_name):
    # Load RHmin data from NetCDF file
    dataset = nc.Dataset(file_path)
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    rhmin = dataset.variables[variable_name][:]
    return lon, lat, rhmin

def plot_rhmin(lon, lat, rhmin_data):
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})

    # Set the extent to match the data domain
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())

    # Add geographical features
    ax.add_feature(cfeature.LAND, edgecolor='black')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')

    # Define RH bins and colors from 0 to 100 in steps of 10
    bounds = list(range(0, 110, 10))
    colors = [
        '#f7fcf0', '#e0f3db', '#ccebc5', '#a8ddb5', '#7bccc4',
        '#4eb3d3', '#2b8cbe', '#0868ac', '#084081', '#081d58'
    ]
    
    cmap = mcolors.ListedColormap(colors)
    norm = mcolors.BoundaryNorm(bounds, len(bounds) - 1)

    rhmin_slice = rhmin_data[0, :, :]
    masked_rhmin = np.ma.masked_invalid(rhmin_slice)

    # Plot RH
    contour = ax.contourf(
        lon, lat, masked_rhmin,
        levels=bounds,
        cmap=cmap, norm=norm, extend='neither', transform=ccrs.PlateCarree()
    )

    # White ocean mask (using inverted land polygons)
    land_shp = shpreader.natural_earth(resolution='10m', category='physical', name='land')
    geometries = list(shpreader.Reader(land_shp).geometries())

    full_extent = Polygon([
        (lon.min(), lat.min()),
        (lon.max(), lat.min()),
        (lon.max(), lat.max()),
        (lon.min(), lat.max()),
        (lon.min(), lat.min())
    ])

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

    tick_labels = [str(b) for b in bounds]
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.05, ticks=bounds)
    cbar.ax.set_yticklabels(tick_labels)
    cbar.set_label('Minimum Relative Humidity (%)')

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlocator = mticker.MultipleLocator(2)
    gl.ylocator = mticker.MultipleLocator(2)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}

    plt.show()

if __name__ == "__main__":
    rhmin_file = fr'C:\AgERA5\RH\RHmin-2m_C3S-glob-agric_AgERA5_{YYYYMMDD}_final-v1.1.area-subset.-20.15.-35.35.nc'
    rhmin_variable = 'RHmin'

    lon, lat, rhmin_data = load_rhmin_data(rhmin_file, rhmin_variable)
    plot_rhmin(lon, lat, rhmin_data)
