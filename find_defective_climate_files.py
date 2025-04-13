# -*- coding: utf-8 -*-
"""
This script will scan all the CSV files and check whetehr the Tmax value for a specified date (e.g. 20241231)
is blank or missing. If this is the case it will add the coordinates to the file AgERA5_GridPoint_Elevations_sample.csv

Created April 2025
@author: steynas@ufs.ac.za
"""

import os
import pandas as pd

# Directory containing the climate CSV files
CLIMATE_DIR = "/home/steynas/AgERA5/ClimDataBase"
ELEVATION_FILE = "/home/steynas/AgERA5/DEM/AgERA5_GridPoint_Elevations.csv"
SAMPLE_OUTPUT = "/home/steynas/AgERA5/DEM/AgERA5_GridPoint_Elevations_sample.csv"
TARGET_DATE = "20241231"

# Load full elevation data
full_elev_df = pd.read_csv(ELEVATION_FILE)

# Store defective coordinate tuples
defective_coords = []

# Loop through all CSV files in the climate data directory
for fname in os.listdir(CLIMATE_DIR):
    if fname.endswith(".csv") and fname.startswith("ClimateTimeSeries_AgERA5"):
        fpath = os.path.join(CLIMATE_DIR, fname)
        try:
            df = pd.read_csv(fpath, dtype=str)

            # Ensure required columns exist
            if "Date" not in df.columns or "Tmax" not in df.columns:
                continue

            # Find 20241231 row
            row = df[df["Date"] == TARGET_DATE]
            if row.empty:
                continue

            tmax_val = row["Tmax"].values[0]
            if pd.isna(tmax_val) or str(tmax_val).strip() == "":
                print(f"Defective file: {fname}")
                parts = fname.split("_")
                lat = float(parts[-2])
                lon = float(parts[-1].replace(".csv", ""))
                defective_coords.append((round(lon, 1), round(lat, 1)))

        except Exception:
            continue

# Create output DataFrame
sample_rows = []
for lon, lat in defective_coords:
    match = full_elev_df[
        (full_elev_df["latitude"].round(1) == lat) &
        (full_elev_df["longitude"].round(1) == lon)
    ]
    if not match.empty:
        row = match.iloc[0]
        sample_rows.append({
            "longitude": row["longitude"],
            "latitude": row["latitude"],
            "elevation": row["elevation"]
        })

sample_df = pd.DataFrame(sample_rows)
sample_df.to_csv(SAMPLE_OUTPUT, index=False)
print(f"[INFO] Sample file written to {SAMPLE_OUTPUT} with {len(sample_df)} entries.")
print(f"[INFO] Total defective files identified: {len(defective_coords)}")
