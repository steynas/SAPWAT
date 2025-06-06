# -*- coding: utf-8 -*-
"""
Script to plot Wind Speed at 2m (U2) from netCDF using FAO56 conversion.
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

def reduce_wind_speed_to_2m_FAO56(wind_speed_10m):
    adjustment_factor = 4.87 / np.log(67.8 * 10 - 5.42)
    return wind_speed_10m * adjustment_factor

def load_ws_data(file_path, variable_name):
    dataset = nc.Dataset(file_path)
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    wind10m = dataset.variables[variable_name][0, :, :]
    wind2m = reduce_wind_speed_to_2m_FAO56(wind10m)
    return lon, lat, wind2m

def plot_wind2m(lon, lat, ws2m_data):
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.LAND, edgecolor='black')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')

    # Define wind speed bins (m/s) and colours
    bounds = [0, 1, 2, 3, 4, 5, 6, 8, 10, 15]
    colours = [
        '#ffffff',  # 0 m/s (white)
        '#e0f3db',  # 0–1
        '#ccebc5',  # 1–2
        '#a8ddb5',  # 2–3
        '#7bccc4',  # 3–4
        '#4eb3d3',  # 4–5
        '#2b8cbe',  # 5–6
        '#0868ac',  # 6–8
        '#084081',  # 8–10
        '#191970'   # 10–15+
    ]

    cmap = mcolors.ListedColormap(colours)
    norm = mcolors.BoundaryNorm(bounds, len(bounds) - 1)

    masked_ws = np.ma.masked_invalid(ws2m_data)

    contour = ax.contourf(
        lon, lat, masked_ws,
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

    tick_labels = ['1', '2', '3', '4', '5', '6', '8', '10', '15']
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.05, ticks=bounds[1:])
    cbar.ax.set_yticklabels(tick_labels)
    cbar.set_label('Wind Speed at 2 m (m/s)')

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlocator = mticker.MultipleLocator(2)
    gl.ylocator = mticker.MultipleLocator(2)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}

    plt.show()

if __name__ == "__main__":
    ws_file = fr'C:\AgERA5\WS\Wind-Speed-10m-Mean_C3S-glob-agric_AgERA5_{YYYYMMDD}_final-v1.1.area-subset.-20.15.-35.35.nc'
    ws_variable = 'Wind_Speed_10m_Mean'

    lon, lat, ws2m = load_ws_data(ws_file, ws_variable)
    plot_wind2m(lon, lat, ws2m)
