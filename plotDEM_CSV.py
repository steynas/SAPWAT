# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 03:39:13 2024

@author: SteynAS
"""

import pandas as pd
import matplotlib.pyplot as plt

# Path to the CSV file in your local directory
file_path = r'C:\AgERA5\DEM\AgERA5_GridPoint_Elevations.csv'

# Load the CSV data
data = pd.read_csv(file_path)

# Extract latitude, longitude, and elevation
lon = data['longitude'].values
lat = data['latitude'].values
elevation = data['elevation'].values

# Plot Elevation
plt.figure(figsize=(10, 6))

# Create a scatter plot
plt.scatter(lon, lat, c=elevation, cmap='terrain', s=10)

# Add color bar and labels
plt.colorbar(label='Elevation (m)')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# Title is omitted as per your request
# plt.title('Elevation Map')

plt.show()
