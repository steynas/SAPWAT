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
instruction_text = None

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
    global selected_lon, selected_lat, popup_text, closest_gridpoint, instruction_text

    # Ignore clicks outside the map (e.g., in the button area)
    if event.inaxes != ax:
        return

    # Remove the instruction message on first click
    if instruction_text:
        instruction_text.remove()
        instruction_text = None

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
    global selected_lon, selected_lat, popup_text, instruction_text

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
    # Namibia
    "Windhoek": (17.08, -22.57),
    "Keetmanshoop": (18.13, -26.57),
    "Mariental": (18.14, -24.63),
    "Ruacana": (14.42, -17.42),
    "Aussenkehr": (17.43, -28.26),
    "Otjiwarongo": (16.65, -20.46),
    "Gobabis": (18.97, -22.45),
    "Lüderitz": (15.16, -26.65),
    # Botswana
    "Gaborone": (25.92, -24.65),
    "Serowe": (26.71, -22.41),
    "Ghanzi": (21.70, -21.70),
    "Francistown": (27.51, -21.17),
    "Tshane": (21.87, -24.02),
    "Maun": (23.42, -20.00),
    "Orapa": (25.39, -21.31),
    # Zimbabwe
    "Bulawayo": (28.58, -20.15),
    "Chiredzi": (31.67, -21.05),
    "Masvingo": (30.82, -20.06),
    "Chipinge": (32.62, -20.20),
    # Mozambique
    "Maputo": (32.59, -25.96),
    "Chokwe": (33.00, -24.53),
    "Xai-Xai": (33.64, -25.04),
    "Mabote": (34.13, -22.04),
    # Eswatini
    "Mbabane": (31.13, -26.32),
    # Lesotho
    "Maseru": (27.48, -29.31),
    # Limpopo
    "Polokwane": (29.47, -23.90),
    "Thohoyandou": (30.48, -22.95),
    "Hoedspruit": (30.96, -24.35),
    "Groblersdal": (29.39, -25.17),
    "Lephalale": (27.74, -23.66),
    "Musina": (30.05, -22.35),
    "Modimolle": (28.41, -24.70),
    # Kwazulu-Natal
    "Pietermaritzburg": (30.38, -29.61),
    "Newcastle": (29.92, -27.75),
    "Ulundi": (31.41, -28.33),
    "Jozini": (32.06, -27.43),
    "Bergville": (29.37, -28.73),
    "Durban": (31.02, -29.85),
    "Richards Bay": (32.08, -28.78),
    "Kokstad": (29.42, -30.55),
    "Port Edward": (30.23, -31.05),
    # Free State
    "Bloemfontein": (26.21, -29.12),
    "Koffiefontein": (25.00, -29.42),
    "Bethlehem": (28.31, -28.23),
    "Kroonstad": (27.66, -27.65),
    "Bothaville": (26.62, -27.39),
    # North-West
    "Mahikeng": (25.65, -25.85),
    "Hartswater": (24.79, -27.78),
    "Potchefstroom": (27.11, -26.72),
    "Vryburg": (24.73, -26.96),
    "Pomfret": (23.22, -25.87),
    # Northern Cape
    "Kimberley": (24.77, -28.74),
    "Kuruman": (23.43, -27.45),
    "Upington": (21.25, -28.45),
    "De Aar": (24.01, -30.65),
    "Vanderkloof": (24.74, -29.99),
    "Brandvlei": (20.47, -30.47),
    "Calvinia": (19.77, -31.47),
    "Douglas": (23.77, -29.06),
    "Kakamas": (20.62, -28.78),
    "Springbok": (17.88, -29.67),
    "Twee Rivieren": (20.61, -26.47),
    # Western Cape
    "Cape Town": (18.42, -33.93),
    "Beaufort West": (22.57, -32.36),
    "Worcester": (19.44, -33.65),
    "Vredenburg": (17.99, -32.91),
    "Bredasdorp": (20.04, -34.53),
    "Riversdale": (21.25, -34.09),
    "Oudtshoorn": (22.20, -33.59),
    "George": (22.45, -33.96),
    "Vredendal": (18.50, -31.67),
    # Eastern Cape
    "Gqeberha": (25.60, -33.96),
    "East London": (27.91, -33.02),
    "Aliwal North": (26.71, -30.69),
    "Bhisho": (27.31, -32.85),
    "Kirkwood": (25.45, -33.39),
    "Cradock": (25.61, -32.16),
    "Queenstown": (26.88, -31.90),
    "Patensie": (24.82, -33.78),
    "Mthatha": (28.78, -31.59),
    # Gauteng
    "Pretoria": (28.22, -25.75),
    # Mpumalanga
    "Mbombela": (31.03, -25.47),
    "Ermelo": (29.99, -26.52)
    }

    for city, coord in cities.items():
        ax.plot(coord[0], coord[1], marker='o', color='red', markersize=3, transform=ccrs.PlateCarree())
        ax.text(coord[0] + 0.2, coord[1] - 0.2, city, transform=ccrs.PlateCarree(), fontsize=6)

    # Add an instruction message at the top
    instruction_text = ax.text(0.5, 1.05, "Click on any land area to select coordinates",
                               transform=ax.transAxes, fontsize=12, ha='center', color='blue')

    # Create a button to lock in the coordinates
    button_ax = fig.add_axes([0.7, 0.05, 0.2, 0.05])  # Bottom-right corner
    button = Button(button_ax, 'Accept selected coordinates')
    button.on_clicked(on_button_clicked)
    button_ax.set_visible(False)  # Hide button initially

    # Connect the click event
    fig.canvas.mpl_connect('button_press_event', lambda event: onclick(event, ax, button))

    # Display the plot
    plt.show()

# Call the function to


# Call the function to plot the clickable map
plot_clickable_map()
