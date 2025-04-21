# -*- coding: utf-8 -*-
"""
Script to plot the difference between ERA5-Land Tmin and AgERA5 Tmin for a specific date.
Difference = ERA5-Land - AgERA5 (°C)
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
from shapely.geometry import Polygon
from shapely.ops import unary_union
import cartopy.io.shapereader as shpreader
import matplotlib.patches as mpatches
from scipy.interpolate import RegularGridInterpolator

# Set target date string
YYYYMMDD = "20231127"

# Define bins and colour map for difference
bounds = [-10, -5, -3, -2, -1, 1, 2, 3, 4, 5, 10]
bin_colors = [
    '#08306b',  # -10 to -5
    '#2171b5',  # -5 to -3
    '#6baed6',  # -3 to -2
    '#c6dbef',  # -2 to -1
    '#ffffff',  # -1 to +1
    '#fee0d2',  # +1 to +2
    '#fcbba1',  # +2 to +3
    '#fb6a4a',  # +3 to +4
    '#cb181d',  # +4 to +5
    '#67000d'   # +5 to +10
]
cmap = mcolors.ListedColormap(bin_colors)
cmap.set_under('#041c3c')
cmap.set_over('#40000b')
norm = mcolors.BoundaryNorm(bounds, ncolors=len(bin_colors))

def load_tmin_data_ag_era5(agera5_path, era5_path):
    # Load AgERA5
    ds_ag = nc.Dataset(agera5_path)
    lon = ds_ag.variables['lon'][:]
    lat = ds_ag.variables['lat'][:]
    tmin_ag = ds_ag.variables['Temperature_Air_2m_Min_24h'][0, :, :] - 273.15
    ds_ag.close()

    # Load ERA5-Land
    ds_era = nc.Dataset(era5_path)
    lon_era = ds_era.variables['longitude'][:]
    lat_era = ds_era.variables['latitude'][:]
    tmin_era = ds_era.variables['t2m'][0, :, :] - 273.15
    ds_era.close()

    # Align ERA5 to AgERA5 grid by interpolating if needed
    if lat_era.shape != lat.shape or lon_era.shape != lon.shape:
        interp_func = RegularGridInterpolator((lat_era[::-1], lon_era), tmin_era[::-1, :], bounds_error=False, fill_value=np.nan)
        lon2d, lat2d = np.meshgrid(lon, lat)
        coords = np.stack([lat2d.ravel(), lon2d.ravel()], axis=-1)
        tmin_era_interp = interp_func(coords).reshape(lat.shape[0], lon.shape[0])
        tmin_era = tmin_era_interp

    tmin_diff = tmin_era - tmin_ag
    return lon, lat, tmin_diff

def plot_difference(lon, lat, tmin_diff):
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())

    lon2d, lat2d = np.meshgrid(lon, lat)
    full_extent = Polygon([
        (lon2d.min(), lat2d.min()), (lon2d.max(), lat2d.min()),
        (lon2d.max(), lat2d.max()), (lon2d.min(), lat2d.max()),
        (lon2d.min(), lat2d.min())
    ])

    land_shp = shpreader.natural_earth(resolution='10m', category='physical', name='land')
    land_geoms = list(shpreader.Reader(land_shp).geometries())
    land_union = unary_union(land_geoms)
    ocean_geom = full_extent.difference(land_union)

    masked_diff = np.ma.masked_invalid(tmin_diff)
    contour = ax.contourf(
        lon, lat, masked_diff,
        levels=bounds, cmap=cmap, norm=norm,
        extend='both', transform=ccrs.PlateCarree(), zorder=2
    )

    if ocean_geom.geom_type == 'Polygon':
        patch = mpatches.Polygon(list(ocean_geom.exterior.coords), facecolor='white', edgecolor='none',
                                 transform=ccrs.PlateCarree(), zorder=3)
        ax.add_patch(patch)
    elif ocean_geom.geom_type == 'MultiPolygon':
        for poly in ocean_geom.geoms:
            patch = mpatches.Polygon(list(poly.exterior.coords), facecolor='white', edgecolor='none',
                                     transform=ccrs.PlateCarree(), zorder=3)
            ax.add_patch(patch)

    ax.add_feature(cfeature.LAND, facecolor='none', edgecolor='black', zorder=4)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linestyle=':', zorder=5)
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none', zorder=5)
    ax.add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=0.8, zorder=6)

    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.05, extend='both', ticks=bounds)
    cbar.set_label('Tmin Difference (ERA5-Land - AgERA5) [°C]')

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlocator = mticker.MultipleLocator(2)
    gl.ylocator = mticker.MultipleLocator(2)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}

    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.2)

    plt.title(f'Tmin Difference (ERA5-Land - AgERA5) on {YYYYMMDD}')
    plt.show()

if __name__ == "__main__":
    agera5_path = fr'C:\AgERA5\Temps\Temperature-Air-2m-Min-24h_C3S-glob-agric_AgERA5_{YYYYMMDD}_final-v1.1.area-subset.-20.15.-35.35.nc'
    era5_path = fr'C:\AgERA5\Temps\Patch\ERA5-Land_Tmin_{YYYYMMDD}_area-subset.-20.15.-35.35.nc'
    lon, lat, tmin_diff = load_tmin_data_ag_era5(agera5_path, era5_path)
    plot_difference(lon, lat, tmin_diff)
