# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 01:06:28 2024

@author: SteynAS
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import Point
from cartopy.io.shapereader import natural_earth, Reader
from matplotlib.widgets import Button
import os
import pandas as pd
from math import radians, cos, sin, sqrt, atan2

# Global variables to store coordinates
selected_lon = None
selected_lat = None
popup_text = None
closest_gridpoint = None

# Path to the AgERA5 DEM data
base_dir = "C:/AgERA5"
dem_dir = os.path.join(base_dir, "DEM")
gridpoint_file = os.path.join(dem_dir, "AgERA5_GridPoint_Elevations.csv")

# Function to load AgERA5 gridpoint data
def load_gridpoints(filepath):
    grid_data = pd.read_csv(filepath)
    return grid_data[['longitude', 'latitude']].dropna()  # Drop any rows with NaN

# Load gridpoint data
gridpoints = load_gridpoints(gridpoint_file)

# Function to calculate the distance between two points using the Haversine formula
def haversine(lon1, lat1, lon2, lat2):
    # Radius of the Earth in kilometers
    R = 6371.0
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))

    distance = R * c
    return distance

# Function to find the closest AgERA5 gridpoint
def find_closest_gridpoint(lon, lat):
    distances = gridpoints.apply(lambda row: haversine(lon, lat, row['longitude'], row['latitude']), axis=1)
    closest_idx = distances.idxmin()
    return gridpoints.iloc[closest_idx]

# Function to check if the clicked coordinates fall over land
def is_point_on_land(lon, lat):
    land_shapefile = natural_earth(resolution='110m', category='physical', name='land')
    reader = Reader(land_shapefile)
    land_geom = list(reader.geometries())  # Get land polygons
    
    # Create a Point object for the clicked coordinates
    point = Point(lon, lat)
    
    # Check if the point is within any of the land polygons
    for polygon in land_geom:
        if polygon.contains(point):
            return True
    return False

# Function to handle the click event on the map
def onclick(event, ax, button):
    global selected_lon, selected_lat, popup_text, closest_gridpoint

    # Ignore clicks outside the map (e.g., in the button area)
    if event.inaxes != ax:
        return

    # Get the clicked coordinates
    lon, lat = event.xdata, event.ydata
    
    # Clear previous popup message
    if popup_text:
        popup_text.remove()
    
    # Check if the clicked point is on land
    if is_point_on_land(lon, lat):
        selected_lon, selected_lat = lon, lat
        popup_text = ax.text(0.5, 1.05, f"Valid coordinates: Longitude: {lon:.2f}, Latitude: {lat:.2f}",
                             transform=ax.transAxes, fontsize=12, ha='center', color='green')
        button.ax.set_visible(True)  # Show the button

        # Find the closest AgERA5 gridpoint
        closest_gridpoint = find_closest_gridpoint(lon, lat)
    else:
        selected_lon, selected_lat = None, None
        popup_text = ax.text(0.5, 1.05, "Invalid coordinates (over ocean). Please try again.",
                             transform=ax.transAxes, fontsize=12, ha='center', color='red')
        button.ax.set_visible(False)  # Hide the button if invalid coordinates

    plt.draw()

# Function to handle the button click (finalize selection)
def on_button_clicked(event):
    global selected_lon, selected_lat, closest_gridpoint
    if selected_lon is not None and selected_lat is not None:
        # Print the selected coordinates
        print(f"Selected Coordinates: Longitude: {selected_lon:.2f}, Latitude: {selected_lat:.2f}")
        
        # Print the closest gridpoint coordinates
        print(f"Closest AgERA5 Gridpoint: Longitude: {closest_gridpoint['longitude']:.2f}, Latitude: {closest_gridpoint['latitude']:.2f}")
        
        plt.close()  # Close the map when the coordinates are confirmed

# Function to plot the map and allow coordinate selection
def plot_clickable_map():
    global selected_lon, selected_lat, popup_text

    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})

    # Set the extent to match the study area (longitudes 15-35, latitudes -35 to -20)
    ax.set_extent([15, 35, -35, -20], crs=ccrs.PlateCarree())

    # Add land and coastlines
    ax.add_feature(cfeature.LAND, edgecolor='black')
    ax.add_feature(cfeature.COASTLINE, linewidth=1)

    # Add provincial boundaries (admin level 1)
    ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=1, edgecolor='black')

    # Add rivers and lakes (using built-in Cartopy features, with thinner rivers)
    ax.add_feature(cfeature.RIVERS, edgecolor='blue', linewidth=0.5)
    ax.add_feature(cfeature.LAKES, edgecolor='blue', facecolor='none', linewidth=1)

    # Plot major cities
    cities = {
        "Cape Town": (18.42, -33.93),
        "Kimberley": (24.77, -28.74),
        "Bloemfontein": (26.21, -29.12),
        "Bhisho": (27.31, -32.85),
        "Pietermaritzburg": (30.38, -29.61),
        "Pretoria": (28.22, -25.75),
        "Polokwane": (29.47, -23.90),
        "Mahikeng": (25.65, -25.85),
        "Mbombela": (31.03, -25.47),
        "Upington": (21.25, -28.45),
        "Beaufort West": (22.57, -32.36),
        "Newcastle": (29.92, -27.75),
        "Calvinia": (19.77, -31.47),
        "Gqeberha": (25.60, -33.96),
        "Vaalharts": (24.85, -27.95),
        "Windhoek": (17.08, -22.57),
        "Keetmanshoop": (18.13, -26.57),
        "Gaborone": (25.92, -24.65),
        "Serowe": (26.71, -22.41),
        "Ghanzi": (21.70, -21.70),
        "Bulawayo": (28.58, -20.15),
        "Chiredzi": (31.67, -21.05),
        "Maputo": (32.59, -25.96),
        "Mbabane": (31.13, -26.32)
    }

    for city, coord in cities.items():
        ax.plot(coord[0], coord[1], marker='o', color='red', markersize=5, transform=ccrs.PlateCarree())
        ax.text(coord[0] + 0.2, coord[1] - 0.2, city, transform=ccrs.PlateCarree(), fontsize=9)

    # Create a button to lock in the coordinates
    button_ax = fig.add_axes([0.7, 0.05, 0.2, 0.05])  # Bottom-right corner
    button = Button(button_ax, 'Accept selected coordinates')
    button.on_clicked(on_button_clicked)
    button_ax.set_visible(False)  # Hide button initially

    # Connect the click event
    fig.canvas.mpl_connect('button_press_event', lambda event: onclick(event, ax, button))

    # Display the plot
    plt.show()

# Call the function to plot the clickable map
plot_clickable_map()
