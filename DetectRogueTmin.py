# -*- coding: utf-8 -*-
"""
Detect rogue Tmin values in AgERA5 files on HPC cluster.
Scans all files matching the AgERA5 Tmin filename pattern.
Flags values below -30 °C and writes to RogueTmin.csv.
Files are sorted chronologically and limited to a start/end date range.
After scanning, a gap analysis report is saved to RogueTmin_GapSummary.csv.
"""

import numpy as np
import netCDF4 as nc
import os
import csv
import glob
import re
from collections import defaultdict
from datetime import datetime

data_dir = "/home/steynas/AgERA5/Temps"
output_dir = "/home/steynas/AgERA5"

def load_temperature_data(file_path, variable_name):
    dataset = nc.Dataset(file_path)
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    temp = dataset.variables[variable_name][:]
    temp_2d = temp[0, :, :] - 273.15
    dataset.close()
    return lon, lat, temp_2d

def write_rogue_tmin_to_csv(lon, lat, temp_data, date_str, output_path, threshold=-30):
    rogue_found = False
    with open(output_path, mode='a', newline='') as file:
        writer = csv.writer(file)
        for i in range(temp_data.shape[0]):
            for j in range(temp_data.shape[1]):
                tmin_c = temp_data[i, j]
                if not np.isnan(tmin_c) and tmin_c <= threshold:
                    writer.writerow([lon[j], lat[i], date_str, round(tmin_c, 2)])
                    rogue_found = True
    return rogue_found

def analyse_gap_lengths(input_csv, output_summary):
    records = defaultdict(list)
    with open(input_csv, mode='r') as file:
        next(file)  # Skip header
        for row in csv.reader(file):
            lon, lat, date_str, _ = row
            key = (float(lat), float(lon))
            records[key].append(date_str)

    with open(output_summary, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Latitude', 'Longitude', 'StartDate', 'EndDate', 'GapLength'])
        for (lat, lon), dates in records.items():
            sorted_dates = sorted(datetime.strptime(d, "%Y%m%d") for d in dates)
            start = sorted_dates[0]
            current = start
            for i in range(1, len(sorted_dates)):
                if (sorted_dates[i] - current).days == 1:
                    current = sorted_dates[i]
                else:
                    writer.writerow([lat, lon, start.strftime("%Y%m%d"), current.strftime("%Y%m%d"), (current - start).days + 1])
                    start = sorted_dates[i]
                    current = start
            writer.writerow([lat, lon, start.strftime("%Y%m%d"), current.strftime("%Y%m%d"), (current - start).days + 1])

def scan_all_files(start_date, end_date):
    def extract_date_from_filename(filename):
        match = re.search(r'_AgERA5_(\d{8})_', filename)
        return match.group(1) if match else "00000000"

    output_csv = os.path.join(output_dir, "RogueTmin.csv")
    output_dirname = os.path.dirname(output_csv)
    if output_dirname:
        os.makedirs(output_dirname, exist_ok=True)

    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Longitude', 'Latitude', 'Date', 'Tmin (°C)'])

    nc_files = sorted(glob.glob(os.path.join(data_dir, 'Temperature-Air-2m-Min-24h_C3S-glob-agric_AgERA5_*.nc')))
    print(f"Files found: {len(nc_files)}")

    for nc_file in nc_files:
        date_str = extract_date_from_filename(nc_file)
        if not (start_date <= date_str <= end_date):
            continue
        lon, lat, temp = load_temperature_data(nc_file, 'Temperature_Air_2m_Min_24h')
        rogue_found = write_rogue_tmin_to_csv(lon, lat, temp, date_str, output_csv)
        print(f"Processed {date_str}")
        if rogue_found:
            print(f"  → Rogue values detected on {date_str}")

    gap_summary_path = os.path.join(output_dir, "RogueTmin_GapSummary.csv")
    analyse_gap_lengths(output_csv, gap_summary_path)
    print(f"Gap analysis written to: {gap_summary_path}")

if __name__ == "__main__":
    start_date = input("Enter start date (YYYYMMDD): ")
    end_date = input("Enter end date (YYYYMMDD): ")
    scan_all_files(start_date, end_date)
