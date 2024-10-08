# -*- coding: utf-8 -*-
"""
This is a script to extract AgERA5 data, elevation data, 
calculated RHmax, RHmin and ETo data and create a long-term
time series for each land grid point over southern Africa and
save it to a .csv file using Multiprocessing Pool and Function Parallelization

Created in October 2024
@author: SteynAS@ufs.ac.za
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
import calendar
import multiprocessing as mp

# Load elevation data from CSV
def load_elevation_data(elevation_file):
    elev_df = pd.read_csv(elevation_file)
    return elev_df

# Generate date range for the specified start and end month (YYYYMM format)
def generate_date_range(start_yyyymm, end_yyyymm):
    start_year = int(start_yyyymm[:4])
    start_month = int(start_yyyymm[4:])
    end_year = int(end_yyyymm[:4])
    end_month = int(end_yyyymm[4:])
    dates = pd.date_range(start=f"{start_year}-{start_month}-01", end=f"{end_year}-{end_month}-{calendar.monthrange(end_year, end_month)[1]}", freq='D')
    return dates

# Load data from .nc files for a specific variable and gridpoint
def extract_variable_data(dates, var_dir, var_name, file_prefix, target_lat, target_lon):
    data_list = []
    for date in dates:
        file_name = os.path.join(var_dir, f"{file_prefix}_C3S-glob-agric_AgERA5_{date.strftime('%Y%m%d')}_final-v1.1.area-subset.-20.15.-35.35.nc")
        if not os.path.exists(file_name):
            data_list.append(np.nan)
            continue
        try:
            ds = xr.open_dataset(file_name)
            lat_idx = np.abs(ds['lat'] - target_lat).argmin()
            lon_idx = np.abs(ds['lon'] - target_lon).argmin()
            data_value = ds[var_name].isel(lat=lat_idx, lon=lon_idx).values[0]
            data_list.append(data_value)
        except:
            data_list.append(np.nan)
    return np.array(data_list)

# Create the .csv file for a single grid point
def process_gridpoint(row, dates, var_dirs, output_dir):
    target_lat = row['latitude']
    target_lon = row['longitude']
    elevation = row['elevation']
    
    if pd.isna(elevation):
        return
    
    # Create output file path for this grid point
    lat_str = f"{target_lat:.1f}"
    lon_str = f"{target_lon:.1f}"
    output_file = os.path.join(output_dir, f"ClimateTimeSeries_AgERA5_{lat_str}_{lon_str}.csv")

    # Check if the CSV file already exists; if so, skip processing this grid point
    if os.path.exists(output_file):
        print(f"Skipping grid point {target_lat}, {target_lon} as CSV already exists.")
        return
    
    # Extract data for each variable
    tmax_data = extract_variable_data(dates, var_dirs["Tmax"], "Temperature_Air_2m_Max_24h", "Temperature-Air-2m-Max-24h", target_lat, target_lon) - 273.15
    tmin_data = extract_variable_data(dates, var_dirs["Tmin"], "Temperature_Air_2m_Min_24h", "Temperature-Air-2m-Min-24h", target_lat, target_lon) - 273.15
    rhmax_data = extract_variable_data(dates, var_dirs["RHmax"], "RHmax", "RHmax-2m", target_lat, target_lon)
    rhmin_data = extract_variable_data(dates, var_dirs["RHmin"], "RHmin", "RHmin-2m", target_lat, target_lon)
    vpmean_data = extract_variable_data(dates, var_dirs["VPmean"], "Vapour_Pressure_Mean", "Vapour-Pressure-Mean", target_lat, target_lon) * 0.1
    wsmean_data = extract_variable_data(dates, var_dirs["WSmean"], "Wind_Speed_10m_Mean", "Wind-Speed-10m-Mean", target_lat, target_lon) * 0.748
    sr_data = extract_variable_data(dates, var_dirs["SR"], "Solar_Radiation_Flux", "Solar-Radiation-Flux", target_lat, target_lon) * 1e-6
    pr_data = extract_variable_data(dates, var_dirs["Pr"], "Precipitation_Flux", "Precipitation-Flux", target_lat, target_lon)
    eto_data = extract_variable_data(dates, var_dirs["ETo"], "Reference_Evapotranspiration_(PM-FAO56)", "ETo", target_lat, target_lon)

    # Create a DataFrame and save to CSV
    df = pd.DataFrame({
        "Longitude": [f"{target_lon:.1f}"] * len(dates),
        "Latitude": [f"{target_lat:.1f}"] * len(dates),
        "Elevation": [f"{elevation:.0f}"] * len(dates),
        "Date": dates.strftime('%Y%m%d'),
        "Tmax": [f"{x:.2f}" if not np.isnan(x) else "" for x in tmax_data],
        "Tmin": [f"{x:.2f}" if not np.isnan(x) else "" for x in tmin_data],
        "RHmax": [f"{x:.2f}" if not np.isnan(x) else "" for x in rhmax_data],
        "RHmin": [f"{x:.2f}" if not np.isnan(x) else "" for x in rhmin_data],
        "VPmean": [f"{x:.2f}" if not np.isnan(x) else "" for x in vpmean_data],
        "WSmean": [f"{x:.2f}" if not np.isnan(x) else "" for x in wsmean_data],
        "SR": [f"{x:.2f}" if not np.isnan(x) else "" for x in sr_data],
        "Pr": [f"{x:.2f}" if not np.isnan(x) else "" for x in pr_data],
        "ETo": [f"{x:.2f}" if not np.isnan(x) else "" for x in eto_data]
    })
    
    os.makedirs(output_dir, exist_ok=True)
    df.to_csv(output_file, index=False)
    print(f"Processed and saved data for grid point at {target_lat}, {target_lon}")

# Main function with multiprocessing
def create_climate_csv_for_gridpoints_parallel(start_yyyymm, end_yyyymm, elevation_file, output_dir):
    elev_df = load_elevation_data(elevation_file)
    dates = generate_date_range(start_yyyymm, end_yyyymm)

    base_dir = "/home/steynas/AgERA5/"
    var_dirs = {
        "Tmax": os.path.join(base_dir, "Temps"),
        "Tmin": os.path.join(base_dir, "Temps"),
        "RHmax": os.path.join(base_dir, "RH"),
        "RHmin": os.path.join(base_dir, "RH"),
        "VPmean": os.path.join(base_dir, "VP"),
        "WSmean": os.path.join(base_dir, "WS"),
        "SR": os.path.join(base_dir, "SR"),
        "Pr": os.path.join(base_dir, "Pr"),
        "ETo": os.path.join(base_dir, "ETo")
    }

    # Multiprocessing setup
    num_cores = mp.cpu_count()  # Use all available cores
    pool = mp.Pool(processes=num_cores)

    # Distribute tasks (grid points) to worker processes
    results = [pool.apply_async(process_gridpoint, args=(row, dates, var_dirs, output_dir)) for index, row in elev_df.iterrows()]

    # Ensure all tasks are completed
    for result in results:
        result.get()  # This will raise any exceptions from the worker processes

    pool.close()
    pool.join()

if __name__ == "__main__":
    #start_yyyymm = input("Enter the start date (YYYYMM): ")
    start_yyyymm = "197901"
    #end_yyyymm = input("Enter the end date (YYYYMM): ")
    end_yyyymm = "202408"
    elevation_file = r"/home/steynas/AgERA5/DEM/AgERA5_GridPoint_Elevations.csv"
    output_dir = r"/home/steynas/AgERA5/ClimDataBase"
    
    create_climate_csv_for_gridpoints_parallel(start_yyyymm, end_yyyymm, elevation_file, output_dir)
