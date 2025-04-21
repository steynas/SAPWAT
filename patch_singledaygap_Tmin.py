# -*- coding: utf-8 -*-
"""
Patch single-day rogue Tmin values using linear interpolation.
Reads RogueTmin_GapSummary.csv and updates ClimateTimeSeries_AgERA5_{lat}_{lon}.csv
Only applies interpolation where GapLength == 1.
"""

import os
import csv
from datetime import datetime, timedelta
import pandas as pd

# Input paths
gap_summary_file = "/home/steynas/AgERA5/RogueTmin_GapSummary.csv"
climate_db_dir = "/home/steynas/AgERA5/ClimateDataBase"

# Function to apply interpolation for a single file
def patch_single_day_gap(lat, lon, gap_date):
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

        if prev_date in df.index and next_date in df.index:
            tmin_interp = (df.loc[prev_date, 'Tmin'] + df.loc[next_date, 'Tmin']) / 2
            df.loc[date, 'Tmin'] = round(tmin_interp, 2)
            df.sort_index(inplace=True)
            df.to_csv(file_path)
            print(f"Patched ({lat}, {lon}) on {gap_date} with {tmin_interp:.2f} °C")
        else:
            print(f"Cannot interpolate ({lat}, {lon}) on {gap_date}: missing neighbours")
    except Exception as e:
        print(f"Error processing ({lat}, {lon}) on {gap_date}: {e}")

# Read gap summary and apply interpolation
def patch_all_single_day_gaps():
    with open(gap_summary_file, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            if row['GapLength'] == '1':
                patch_single_day_gap(float(row['Latitude']), float(row['Longitude']), row['StartDate'])

if __name__ == "__main__":
    patch_all_single_day_gaps()
