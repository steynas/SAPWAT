# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 12:56:36 2024

@author: SteynAS
"""

import os
import calendar
import datetime
import numpy as np
import pandas as pd
import xarray as xr

def calculate_julian_day(year, month, day):
    year = int(year)
    month = int(month)
    day = int(day)
    date = datetime.date(year, month, day)
    return date.timetuple().tm_yday

def read_z(z_file, min_lon, max_lon, min_lat, max_lat, lons_vp, lats_vp):
    # Load the elevation data from the CSV file
    z_df = pd.read_csv(z_file)
    # Extract longitude, latitude, and elevation columns
    z_lons = z_df['longitude'].values.round(4)  # Round to 4 decimal places
    z_lats = z_df['latitude'].values.round(4)   # Round to 4 decimal places
    elevations = z_df['elevation'].values
    # Create an empty elevation grid with the shape of lats_vp and lons_vp (lat, lon order)
    z = np.full((len(lats_vp), len(lons_vp)), np.nan, dtype=np.float64)
    # Create a dictionary to map (lat, lon) pairs to elevation values (we map lat first)
    elevation_map = {(lat, lon): elev for lat, lon, elev in zip(z_lats, z_lons, elevations)}
    lons_vp = lons_vp.round(4)  # Round to 4 decimal places
    lats_vp = lats_vp.round(4)  # Round to 4 decimal places
    # Fill the grid based on lats_vp and lons_vp
    for i, lat in enumerate(lats_vp):
        for j, lon in enumerate(lons_vp):
            z[i, j] = elevation_map.get((lat, lon), np.nan)
    return z

def calculate_rhmin_rhmax(year, month, day, rh_dir, min_lat, max_lat, min_lon, max_lon):
    # Variables in RH file: ['Relative_Humidity_2m_06h', 'time', 'lon', 'lat']
    target_date = f"{year}{month.zfill(2)}{str(day).zfill(2)}"
    RH_min_values = None
    RH_max_values = None
    time_to_variable = {
        '06:00': 'Relative_Humidity_2m_06h',
        '09:00': 'Relative_Humidity_2m_09h',
        '12:00': 'Relative_Humidity_2m_12h',
        '15:00': 'Relative_Humidity_2m_15h',
        '18:00': 'Relative_Humidity_2m_18h'
    }
    for filename in sorted(os.listdir(rh_dir)):
        if filename.endswith(".nc") and f"_{target_date}_" in filename:
            file_path = os.path.join(rh_dir, filename)
            ds = xr.open_dataset(file_path)
            for time, variable_name in time_to_variable.items():
                if variable_name in ds:
                    lats = ds['lat'].values
                    lons = ds['lon'].values
                    lat_indices = np.where((lats >= min_lat) & (lats <= max_lat))[0]
                    lon_indices = np.where((lons >= min_lon) & (lons <= max_lon))[0]
                    RH_data = ds[variable_name].values[0, lat_indices.min():lat_indices.max()+1, lon_indices.min():lon_indices.max()+1]
                    if RH_min_values is None:
                        RH_min_values = np.full(RH_data.shape, np.inf)
                        RH_max_values = np.full(RH_data.shape, -np.inf)
                    RH_min_values = np.fmin(RH_min_values, np.nan_to_num(RH_data, nan=np.inf))
                    RH_max_values = np.fmax(RH_max_values, np.nan_to_num(RH_data, nan=-np.inf))
                    break
    if RH_min_values is None or RH_max_values is None:
        print(f"No valid RH data found for {year}-{month.zfill(2)}-{str(day).zfill(2)}. Skipping this day.")
        return None, None, None, None, None, None, None
    RH_min_values[RH_min_values == np.inf] = np.nan
    RH_max_values[RH_max_values == -np.inf] = np.nan
    return RH_min_values, RH_max_values, lats, lons, lat_indices, lon_indices, ds

def read_vp(year, month, day, vp_dir, min_lat, max_lat, min_lon, max_lon):
    # Variables in VP file: ['Vapour_Pressure_Mean', 'time', 'lon', 'lat']
    target_date = f"{year}{month.zfill(2)}{str(day).zfill(2)}"
    vp_file = os.path.join(vp_dir, f"Vapour-Pressure-Mean_C3S-glob-agric_AgERA5_{target_date}_final-v1.1.area-subset.-20.15.-35.35.nc")
    if not os.path.exists(vp_file):
        print(f"Missing VP file for {year}-{month.zfill(2)}-{str(day).zfill(2)}. Skipping this day.")
        return None, None, None, None, None
    vp_ds = xr.open_dataset(vp_file)
    # Extract lat, lon, and time variables
    lats = vp_ds['lat'].values
    lons = vp_ds['lon'].values
    time = vp_ds['time'].values  
    # Get indices for the bounding box
    lat_indices = np.where((lats >= min_lat) & (lats <= max_lat))[0]
    lon_indices = np.where((lons >= min_lon) & (lons <= max_lon))[0]
    # Slice the VP data based on lat/lon indices
    VP = vp_ds['Vapour_Pressure_Mean'].values[:, lat_indices.min():lat_indices.max()+1, lon_indices.min():lon_indices.max()+1]
    return VP, time, lats[lat_indices], lons[lon_indices], lat_indices, lon_indices

def read_tmin_tmax(year, month, day, temp_dir, min_lat, max_lat, min_lon, max_lon):
    # Variables in Tmin file: ['Temperature_Air_2m_Min_24h', 'time', 'lon', 'lat']
    # Variables in Tmax file: ['Temperature_Air_2m_Max_24h', 'time', 'lon', 'lat']
    target_date = f"{year}{month.zfill(2)}{str(day).zfill(2)}"
    Tmin_file = os.path.join(temp_dir, f"Temperature-Air-2m-Min-24h_C3S-glob-agric_AgERA5_{target_date}_final-v1.1.area-subset.-20.15.-35.35.nc")
    Tmax_file = os.path.join(temp_dir, f"Temperature-Air-2m-Max-24h_C3S-glob-agric_AgERA5_{target_date}_final-v1.1.area-subset.-20.15.-35.35.nc")
    if not os.path.exists(Tmin_file) or not os.path.exists(Tmax_file):
        print(f"Missing Tmin or Tmax files for {year}-{month.zfill(2)}-{str(day).zfill(2)}. Skipping this day.")
        return None, None
    # Open both datasets
    Tmin_ds = xr.open_dataset(Tmin_file)
    Tmax_ds = xr.open_dataset(Tmax_file)
    # Extract lat, lon, and time variables
    lats_temp = Tmin_ds['lat'].values
    lons_temp = Tmin_ds['lon'].values
    time_temp = Tmin_ds['time'].values
    # Get indices for the bounding box
    lon_indices_temp = np.where((lons_temp >= min_lon) & (lons_temp <= max_lon))[0]
    lat_indices_temp = np.where((lats_temp >= min_lat) & (lats_temp <= max_lat))[0]
    # Slice the Tmin, Tmax data based on lat/lon indices
    Tmin = Tmin_ds['Temperature_Air_2m_Min_24h'].values[:, lat_indices_temp.min():lat_indices_temp.max()+1, lon_indices_temp.min():lon_indices_temp.max()+1]
    Tmax = Tmax_ds['Temperature_Air_2m_Max_24h'].values[:, lat_indices_temp.min():lat_indices_temp.max()+1, lon_indices_temp.min():lon_indices_temp.max()+1]
    return Tmin, Tmax, time_temp, lats_temp[lat_indices_temp], lons_temp[lon_indices_temp], lat_indices_temp, lon_indices_temp

def read_sr(year, month, day, sr_dir, min_lat, max_lat, min_lon, max_lon):
    # Variables in SR file: ['Solar_Radiation_Flux', 'time', 'lon', 'lat']
    target_date = f"{year}{month.zfill(2)}{str(day).zfill(2)}"
    sr_file = os.path.join(sr_dir, f"Solar-Radiation-Flux_C3S-glob-agric_AgERA5_{target_date}_final-v1.1.area-subset.-20.15.-35.35.nc")
    if not os.path.exists(sr_file):
        print(f"Missing SR file for {year}-{month.zfill(2)}-{str(day).zfill(2)}. Skipping this day.")
        return None, None, None, None, None
    sr_ds = xr.open_dataset(sr_file)
    # Extract lat, lon, and time variables
    lats_sr = sr_ds['lat'].values
    lons_sr = sr_ds['lon'].values
    time_sr = sr_ds['time'].values  
    # Get indices for the bounding box
    lon_indices_sr = np.where((lons_sr >= min_lon) & (lons_sr <= max_lon))[0]
    lat_indices_sr = np.where((lats_sr >= min_lat) & (lats_sr <= max_lat))[0]
    # Slice the SR data based on lat/lon indices
    SR = sr_ds['Solar_Radiation_Flux'].values[:, lat_indices_sr.min():lat_indices_sr.max()+1, lon_indices_sr.min():lon_indices_sr.max()+1]
    return SR, time_sr, lats_sr[lat_indices_sr], lons_sr[lon_indices_sr], lat_indices_sr, lon_indices_sr

def read_ws(year, month, day, ws_dir, min_lat, max_lat, min_lon, max_lon):
    # Variables in WS file: ['Wind_Speed_10m_Mean', 'time', 'lon', 'lat']
    target_date = f"{year}{month.zfill(2)}{str(day).zfill(2)}"
    ws_file = os.path.join(ws_dir, f"Wind-Speed-10m-Mean_C3S-glob-agric_AgERA5_{target_date}_final-v1.1.area-subset.-20.15.-35.35.nc")
    if not os.path.exists(ws_file):
        print(f"Missing WS file for {year}-{month.zfill(2)}-{str(day).zfill(2)}. Skipping this day.")
        return None
    ws_ds = xr.open_dataset(ws_file)
    # Extract lat, lon, and time variables
    lats_ws = ws_ds['lat'].values
    lons_ws = ws_ds['lon'].values
    time_ws = ws_ds['time'].values
    # Get indices for the bounding box
    lon_indices_ws = np.where((lons_ws >= min_lon) & (lons_ws <= max_lon))[0]
    lat_indices_ws = np.where((lats_ws >= min_lat) & (lats_ws <= max_lat))[0]
    # Slice the WS data based on lat/lon indices
    WS = ws_ds['Wind_Speed_10m_Mean'].values[:, lat_indices_ws.min():lat_indices_ws.max()+1, lon_indices_ws.min():lon_indices_ws.max()+1]
    return WS, time_ws, lats_ws[lat_indices_ws], lons_ws[lon_indices_ws], lat_indices_ws, lon_indices_ws

def read_pr(year, month, day, pr_dir, min_lat, max_lat, min_lon, max_lon):
    # Variables in Pr file: ['Precipitation_Flux', 'time', 'lon', 'lat']
    target_date = f"{year}{month.zfill(2)}{str(day).zfill(2)}"
    pr_file = os.path.join(pr_dir, f"Precipitation-Flux_C3S-glob-agric_AgERA5_{target_date}_final-v1.1.area-subset.-20.15.-35.35.nc")
    if not os.path.exists(pr_file):
        print(f"Missing Pr file for {year}-{month.zfill(2)}-{str(day).zfill(2)}. Skipping this day.")
        return None
    pr_ds = xr.open_dataset(pr_file)
    # Extract lat, lon, and time variables
    lats_pr = pr_ds['lat'].values
    lons_pr = pr_ds['lon'].values
    time_pr = pr_ds['time'].values
    # Get indices for the bounding box
    lon_indices_pr = np.where((lons_pr >= min_lon) & (lons_pr <= max_lon))[0]
    lat_indices_pr = np.where((lats_pr >= min_lat) & (lats_pr <= max_lat))[0]
    # Slice the WS data based on lat/lon indices
    Pr = pr_ds['Precipitation_Flux'].values[:, lat_indices_ws.min():lat_indices_ws.max()+1, lon_indices_ws.min():lon_indices_ws.max()+1]
    return Pr, time_pr, lats_pr[lat_indices_pr], lons_pr[lon_indices_pr], lat_indices_pr, lon_indices_pr

    # Expand z to match the (time, lat, lon) dimensions
    z_expanded = z[np.newaxis, :, :]  # Expanding z to have a time dimension
    # Define constants
    G = 0.0                # Soil heat flux density (MJ/m²/day)
    Gsc = 0.0820          # Solar constant (MJ/m²/min)
    sigma = 4.903e-9      # Stefan-Boltzmann constant (MJ/m²/K⁴/day)
    cp = 1.013e-3         # Specific heat capacity of air (MJ/kg/°C)
    epsilon = 0.622       # Ratio of molecular weight of water vapour/dry air
    Lv = 2.45             # Latent heat of vaporization (MJ/kg)
    alpha = 0.23          # Albedo for hypothetical grass reference crop
    # Unit conversions
    phi = np.radians(lats)  # Latitude from ° to rad
    TminC = Tmin - 273.15   # Convert from K to °C
    TmaxC = Tmax - 273.15   # Convert from K to °C
    Rs = SR * 1e-6          # Convert from J/m²/day to MJ/m²/day
    U2 = WS * 0.748         # Convert wind speed from 10 m to 2 m
    ea = VP * 0.1           # Convert from hPa to kPa
    # Mean temperature
    Tmean = (Tmin + Tmax) / 2.0      # (K)
    TmeanC = (TminC + TmaxC) / 2.0   # (°C)
    # Atmospheric pressure (kPa)
    P = 101.3 * ((Tmean - 0.0065 * z_expanded) / Tmean) ** 5.26   # z_expanded is 3D now
    # Psychrometric constant (kPa/°C)
    gamma = (cp * P) / (epsilon * Lv)
    # Saturation vapour pressure (kPa)
    es_Tmin = 0.6108 * np.exp((17.27 * TminC) / (TminC + 237.3))
    es_Tmax = 0.6108 * np.exp((17.27 * TmaxC) / (TmaxC + 237.3))
    es = (es_Tmin + es_Tmax) / 2
    # Slope of the vapour pressure curve (kPa/°C)
    delta = (4098 * es) / ((TmeanC + 237.3) ** 2)
    # Julian day
    J = calculate_julian_day(year, month, day)
    # Extraterrestrial radiation (Ra)
    dr = 1 + 0.033 * np.cos(2 * np.pi * J / 365)
    declination = 0.409 * np.sin((2 * np.pi * J / 365) - 1.39)
    omega = np.arccos(-np.tan(phi) * np.tan(declination))
    Ra = (24 * 60 / np.pi) * Gsc * dr * (omega * np.sin(phi) * 
         np.sin(declination) + np.cos(phi) * np.cos(declination) *
         np.sin(omega))
    # Clear-sky solar radiation (Rso)
    Rso = (0.75 + 2e-5 * z_expanded) * Ra
    # Net shortwave radiation (Rns)
    Rns = (1 - alpha) * Rs
    # Net longwave radiation (Rnl)
    Rnl = sigma * ((Tmax**4 + Tmin**4) / 2) * (0.34 - 0.14 * 
          np.sqrt(ea)) * (1.35 * Rs / Rso - 0.35)
    # Net radiation (Rn)
    Rn = Rns - Rnl
    # Reference evapotranspiration (mm/day) according to FAO56 Penman-Monteith equation
    energy_component = 0.408 * delta * (Rn - G)
    aerodynamic_component = gamma * (900 / (TmeanC + 273)) * U2 * (es - ea)
    denominator = delta + gamma * (1 + 0.34 * U2)
    ETo = (energy_component + aerodynamic_component) / denominator
    # Ensure ETo is a float32 array and add a time dimension (1, lat, lon)
    return np.expand_dims(ETo.astype(np.float32), axis=0)
def calculate_eto(year, month, day, min_lat, max_lat, min_lon, max_lon, lats, lons, z, Tmin, Tmax, VP, SR, WS):
    # Expand z to match the (time, lat, lon) dimensions
    z_expanded = z[np.newaxis, :, :]  # Expanding z to have a time dimension
    # Define constants
    G = 0                 # Soil heat flux density (MJ/m²/day)
    Gsc = 0.0820          # Solar constant (MJ/m²/min)
    sigma = 4.903e-9      # Stefan-Boltzmann constant (MJ/m²/K⁴/day)
    cp = 1.013e-3         # Specific heat capacity of air (MJ/kg/°C)
    epsilon = 0.622       # Ratio of molecular weight of water vapour/dry air
    Lv = 2.45             # Latent heat of vaporization (MJ/kg)
    alpha = 0.23          # Albedo for hypothetical grass reference crop
    # Create 2D meshgrid of latitudes and longitudes
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    # Convert latitudes to radians
    phi = np.radians(lat_grid)  # Latitude from degrees to radians, shape (lat, lon)
    # Unit conversions for other variables
    TminC = Tmin - 273.15   # Convert from K to °C
    TmaxC = Tmax - 273.15   # Convert from K to °C
    Rs = SR * 1e-6          # Convert from J/m²/day to MJ/m²/day
    U2 = WS * 0.748         # Convert wind speed from 10 m to 2 m
    ea = VP * 0.1           # Convert from hPa to kPa
    # Mean temperature (K and °C)
    Tmean = (Tmin + Tmax) / 2      # (K)
    TmeanC = (TminC + TmaxC) / 2   # (°C)
    # Atmospheric pressure (kPa), z_expanded is already (1, lat, lon)
    P = 101.3 * ((Tmean - 0.0065 * z_expanded) / Tmean) ** 5.26
    # Psychrometric constant (kPa/°C)
    gamma = (cp * P) / (epsilon * Lv)
    # Saturation vapour pressure (kPa)
    es_Tmin = 0.6108 * np.exp((17.27 * TminC) / (TminC + 237.3))
    es_Tmax = 0.6108 * np.exp((17.27 * TmaxC) / (TmaxC + 237.3))
    es = (es_Tmin + es_Tmax) / 2
    # Slope of the vapour pressure curve (kPa/°C)
    delta = (4098 * es) / ((TmeanC + 237.3) ** 2)
    # Julian day
    J = calculate_julian_day(year, month, day)
    # Extraterrestrial radiation (Ra) calculation
    dr = 1 + 0.033 * np.cos(2 * np.pi * J / 365)
    declination = 0.409 * np.sin((2 * np.pi * J / 365) - 1.39)
    omega = np.arccos(-np.tan(phi) * np.tan(declination))
    # Compute Ra over the 2D grid of latitudes and longitudes
    Ra = (24 * 60 / np.pi) * Gsc * dr * (
        omega * np.sin(phi) * np.sin(declination)
        + np.cos(phi) * np.cos(declination) * np.sin(omega)
    )
    # Expand Ra to (1, lat, lon) to match the (time, lat, lon) shape
    Ra = np.expand_dims(Ra, axis=0)
    # Clear-sky solar radiation (Rso)
    Rso = (0.75 + 2e-5 * z_expanded) * Ra
    # Net shortwave radiation (Rns)
    Rns = (1 - alpha) * Rs
    # Net longwave radiation (Rnl)
    Rnl = sigma * ((Tmax**4 + Tmin**4) / 2) * (0.34 - 0.14 * np.sqrt(ea)) * (1.35 * Rs / Rso - 0.35)
    # Net radiation (Rn)
    Rn = Rns - Rnl
    # Reference evapotranspiration (mm/day) according to FAO56 Penman-Monteith equation
    energy_component = 0.408 * delta * (Rn - G)
    aerodynamic_component = gamma * (900 / (TmeanC + 273)) * U2 * (es - ea)
    denominator = delta + gamma * (1 + 0.34 * U2)
    ETo = (energy_component + aerodynamic_component) / denominator
    # Ensure ETo is a float32 array with shape (1, lat, lon)
    return ETo.astype(np.float32)

#######################################
## Main script execution starts here ##
#######################################
year = input("Enter the year (e.g., 1979): ")
month = input("Enter the month (01-12): ")
days_in_month = calendar.monthrange(int(year), int(month))[1]

base_dir = "C:/AgERA5"
dem_dir = os.path.join(base_dir, "DEM")
rh_dir = os.path.join(base_dir, "RH")
vp_dir = os.path.join(base_dir, "VP")
temp_dir = os.path.join(base_dir, "Temps")
ws_dir = os.path.join(base_dir, "WS")
sr_dir = os.path.join(base_dir, "SR")
pr_dir = os.path.join(base_dir, "Pr")
eto_dir = os.path.join(base_dir, "ETo")
os.makedirs(eto_dir, exist_ok=True)

# Defined boundaries of geographical area of validity 
min_lat, max_lat = -35, -20
min_lon, max_lon = 15, 35

# Provide elevation data file path and initialize grid
z_file = os.path.join(dem_dir, "AgERA5_GridPoint_Elevations.csv")
z = None

for day in range(1, days_in_month + 1):
    
    # Calculate RHmin and RHmax and save these to individual files
    # RHmin and RHmax are not used in ETo calculation since VP is available
    # However, these values are required for the database
    # Calculate RHmin and RHmax
    RHmin, RHmax, lats, lons, lat_indices, lon_indices, ds = calculate_rhmin_rhmax(
        year, month, day, rh_dir, min_lat, max_lat, min_lon, max_lon)
    if RHmin is None or RHmax is None:
        continue
    # Expand RHmin and RHmax to include the time dimension (still in lat-lon order)
    RHmin = np.expand_dims(RHmin, axis=0)  # Shape becomes (1, 150, 200)
    RHmax = np.expand_dims(RHmax, axis=0)  # Shape becomes (1, 150, 200)
    # Create the time coordinate for this day
    time = np.array([f"{year}-{month.zfill(2)}-{str(day).zfill(2)}T00:00:00"], dtype="datetime64[ns]")
    # Create the Dataset for RHmin with (time, lat, lon)
    output_ds_min = xr.Dataset(
        {"RHmin": (["time", "lat", "lon"], RHmin)},
        coords={
            "time": time,
            "lon": lons[lon_indices],  # Longitude
            "lat": lats[lat_indices]   # Latitude
        },
        attrs=ds.attrs
    )
    # Create the Dataset for RHmax with (time, lat, lon)
    output_ds_max = xr.Dataset(
        {"RHmax": (["time", "lat", "lon"], RHmax)},
        coords={
            "time": time,
            "lon": lons[lon_indices],  # Longitude
            "lat": lats[lat_indices]   # Latitude
        },
        attrs=ds.attrs
    )
    # Save RHmin and RHmax datasets to NetCDF
    target_date = f"{year}{month.zfill(2)}{str(day).zfill(2)}"
    rhmin_filename = os.path.join(rh_dir, f"RHmin-2m_C3S-glob-agric_AgERA5_{target_date}_final-v1.1.area-subset.-20.15.-35.35.nc")
    rhmax_filename = os.path.join(rh_dir, f"RHmax-2m_C3S-glob-agric_AgERA5_{target_date}_final-v1.1.area-subset.-20.15.-35.35.nc")
    output_ds_min.to_netcdf(rhmin_filename, engine="netcdf4")
    output_ds_max.to_netcdf(rhmax_filename, engine="netcdf4")
    
    # Read Vapour Pressure (VP)
    VP, time_vp, lats_vp, lons_vp, lat_indices_vp, lon_indices_vp = read_vp(
        year, month, day, vp_dir, min_lat, max_lat, min_lon, max_lon)
    if VP is None:
        continue

    # Read Tmin and Tmax and convert to degrees Celsius
    Tmin, Tmax, time_temp, lats_temp, lons_temp, lat_indices_temp, lon_indices_temp = read_tmin_tmax(
        year, month, day, temp_dir, min_lat, max_lat, min_lon, max_lon)
    if Tmin is None or Tmax is None:
        continue
    
    # Read Solar Radiation (SR)
    SR, time_sr, lats_sr, lons_sr, lat_indices_sr, lon_indices_sr = read_sr(year, month, day, sr_dir, min_lat, max_lat, min_lon, max_lon)
    if SR is None:
        continue

    # Read mean 10m Wind Speed (WS)
    WS, time_ws, lats_ws, lons_ws, lat_indices_ws, lon_indices_ws = read_ws(year, month, day, ws_dir, min_lat, max_lat, min_lon, max_lon)
    if WS is None:
        continue

    # Read Precipitation (Pr)
    Pr, time_pr, lats_pr, lons_pr, lat_indices_pr, lon_indices_pr = read_pr(year, month, day, pr_dir, min_lat, max_lat, min_lon, max_lon)
    if Pr is None:
        continue
    
    # Extract elevation data for each gridpoint (once only)
    if z is None:  # Check if z has already been calculated
        z = read_z(z_file, min_lon, max_lon, min_lat, max_lat, lons_vp, lats_vp)
        
    # Calculate ETo
    ETo = calculate_eto(year, month, day, min_lat, max_lat, min_lon, max_lon, lats, lons, z, Tmin, Tmax, VP, SR, WS)

    # Create the time coordinate for this day
    time = np.array([f"{year}-{month.zfill(2)}-{str(day).zfill(2)}T00:00:00"], dtype="datetime64[ns]")

    # Create the Dataset for ETo with (time, lat, lon) dimensions
    eto_dataset = xr.Dataset(
        {"Reference_Evapotranspiration_(PM-FAO56)": (["time", "lat", "lon"], ETo)},
        coords={
            "time": time,
            "lon": lons_vp,  # Longitude array
            "lat": lats_vp   # Latitude array
        },
        attrs={"description": "Reference Evapotranspiration (PM-FAO56)"}
    )

    # Define the output file path for saving ETo
    eto_file = os.path.join(eto_dir, f"ETo_C3S-glob-agric_AgERA5_{year}{month.zfill(2)}{str(day).zfill(2)}_final-v1.1.area-subset.-20.15.-35.35.nc")
    
    # Save the ETo dataset to NetCDF
    eto_dataset.to_netcdf(eto_file, engine="netcdf4")
    
    # Message that shows progress
    print(f"Processed {year}-{month.zfill(2)}-{str(day).zfill(2)} successfully.")