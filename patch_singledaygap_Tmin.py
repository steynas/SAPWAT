# -*- coding: utf-8 -*-
"""
Script to patch single-day rogue Tmin-values using linear interpolation.
This script reads RogueTmin_GapSummary.csv and updates ClimateTimeSeries_AgERA5_{lat}_{lon}.csv in ClimDataBase.
Only applies interpolation where GapLength == 1.
Outputs a patch log file: PatchedTmin_Log.csv

Created in April 2025
@author: SteynAS@ufs.ac.za
"""

import os
import csv
from datetime import datetime, timedelta
import pandas as pd

# Input paths
gap_summary_file = "/home/steynas/AgERA5/RogueTmin_GapSummary.csv"
climate_db_dir = "/home/steynas/AgERA5/ClimDataBase"
log_file = "/home/steynas/AgERA5/PatchedTmin_Log.csv"

# Initialise log
def init_log():
    with open(log_file, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Latitude', 'Longitude', 'Date', 'Interpolated Tmin (degC)'])

# Function to apply interpolation for a single file
def patch_single_day_gap(lat, lon, gap_date):
    lat = round(lat, 1)
    lon = round(lon, 1)
    lat_str = f"{lat:.1f}"
    lon_str = f"{lon:.1f}"
    file_name = f"ClimateTimeSeries_AgERA5_{lat_str}_{lon_str}.csv"
    file_path = os.path.join(climate_db_dir, file_name)

    if not os.path.isfile(file_path):
        print(f"Missing file for ({lat}, {lon}): {file_name}")
        return

    df = pd.read_csv(file_path, parse_dates=['Date'])
    df.set_index('Date', inplace=True)

    try:
        date = datetime.strptime(gap_date, "%Y%m%d")
        prev_date = date - timedelta(days=1)
        next_date = date + timedelta(days=1)

        if date not in df.index:
            print(f"Date {gap_date} not found in index for ({lat}, {lon})")
            return

        if prev_date in df.index and next_date in df.index:
            tmin_interp = (df.loc[prev_date, 'Tmin'] + df.loc[next_date, 'Tmin']) / 2
            df.loc[date, 'Tmin'] = round(tmin_interp, 2)
            df.sort_index(inplace=True)
            df_reset = df.reset_index()
            df_reset['Date'] = df_reset['Date'].dt.strftime('%Y%m%d')
            column_order = ['Longitude', 'Latitude', 'Elevation', 'Date', 'Tmax', 'Tmin', 'RHmax', 'RHmin', 'VPmean', 'WSmean', 'SR', 'Pr', 'ETo']
            df_reset = df_reset[column_order]  # Ensure column order
            df_reset.to_csv(file_path, index=False)
            print(f"Patched ({lat}, {lon}) on {gap_date} with {tmin_interp:.2f} °C")

            # Log the patch
            with open(log_file, mode='a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([lat, lon, gap_date, round(tmin_interp, 2)])
        else:
            print(f"Cannot interpolate ({lat}, {lon}) on {gap_date}: missing neighbours")
    except Exception as e:
        print(f"Error processing ({lat}, {lon}) on {gap_date}: {e}")

# Read gap summary and apply interpolation
def patch_all_single_day_gaps():
    init_log()
    patched_files = set()
    with open(gap_summary_file, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            if row['GapLength'] == '1':
                lat = float(row['Latitude'])
                lon = float(row['Longitude'])
                gap_date = row['StartDate']
                lat_rounded = round(lat, 1)
                lon_rounded = round(lon, 1)
                patch_single_day_gap(lat, lon, gap_date)
                patched_files.add((lat_rounded, lon_rounded))
    print(f"\n✅ Patched Tmin in {len(patched_files)} files using linear interpolation.")

if __name__ == "__main__":
    patch_all_single_day_gaps()
