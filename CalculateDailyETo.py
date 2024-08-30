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

def get_elevation(z, lats, lons, elevation_data):
    for i, lat in enumerate(lats):
        for j, lon in enumerate(lons):
            # Round latitude and longitude to match the DEM precision
            rounded_lat = round(lat, 4)
            rounded_lon = round(lon, 4)
            # Match with elevation data
            matching_rows = elevation_data[
                (elevation_data['latitude'].round(4) == rounded_lat) & 
                (elevation_data['longitude'].round(4) == rounded_lon)
            ]
            if not matching_rows.empty:
                z_value = matching_rows['elevation'].values[0]
                z[i, j] = z_value  # Assign elevation value
    return z

def calculate_rhmin_rhmax(year, month, day, rh_dir, min_lat, max_lat, min_lon, max_lon):
    # Variables in RH file: ['Relative_Humidity_2m_06h', 'time', 'lon', 'lat']
    # Variables in RH file: ['Relative_Humidity_2m_09h', 'time', 'lon', 'lat']
    # Variables in RH file: ['Relative_Humidity_2m_12h', 'time', 'lon', 'lat']
    # Variables in RH file: ['Relative_Humidity_2m_15h', 'time', 'lon', 'lat']
    # Variables in RH file: ['Relative_Humidity_2m_18h', 'time', 'lon', 'lat']
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
        return None, None
    RH_min_values[RH_min_values == np.inf] = np.nan
    RH_max_values[RH_max_values == -np.inf] = np.nan
    return RH_min_values, RH_max_values, lats, lons, lat_indices, lon_indices, ds

def read_vp(year, month, day, vp_dir, min_lat, max_lat, min_lon, max_lon):
    # Variables in VP file: ['Vapour_Pressure_Mean', 'time', 'lon', 'lat']
    target_date = f"{year}{month.zfill(2)}{str(day).zfill(2)}"
    vp_file = os.path.join(vp_dir, f"Vapour-Pressure-Mean_C3S-glob-agric_AgERA5_{target_date}_final-v1.1.area-subset.-20.15.-35.35.nc")
    if not os.path.exists(vp_file):
        print(f"Missing VP file for {year}-{month.zfill(2)}-{str(day).zfill(2)}. Skipping this day.")
        return None
    vp_ds = xr.open_dataset(vp_file)
    lats = vp_ds['lat'].values
    lons = vp_ds['lon'].values
    lat_indices = np.where((lats >= min_lat) & (lats <= max_lat))[0]
    lon_indices = np.where((lons >= min_lon) & (lons <= max_lon))[0]
    VP = vp_ds['Vapour_Pressure_Mean'].values[0, lat_indices.min():lat_indices.max()+1, lon_indices.min():lon_indices.max()+1]
    return VP, lats, lons, lat_indices, lon_indices

def read_tmin_tmax(year, month, day, temp_dir, min_lat, max_lat, min_lon, max_lon):
    # Variables in Tmin file: ['Temperature_Air_2m_Min_24h', 'time', 'lon', 'lat']
    # Variables in Tmax file: ['Temperature_Air_2m_Max_24h', 'time', 'lon', 'lat']
    target_date = f"{year}{month.zfill(2)}{str(day).zfill(2)}"
    Tmin_file = os.path.join(temp_dir, f"Temperature-Air-2m-Min-24h_C3S-glob-agric_AgERA5_{target_date}_final-v1.1.area-subset.-20.15.-35.35.nc")
    Tmax_file = os.path.join(temp_dir, f"Temperature-Air-2m-Max-24h_C3S-glob-agric_AgERA5_{target_date}_final-v1.1.area-subset.-20.15.-35.35.nc")
    if not os.path.exists(Tmin_file) or not os.path.exists(Tmax_file):
        print(f"Missing Tmin or Tmax files for {year}-{month.zfill(2)}-{str(day).zfill(2)}. Skipping this day.")
        return None, None
    Tmin_ds = xr.open_dataset(Tmin_file)
    Tmax_ds = xr.open_dataset(Tmax_file)
    lats = Tmin_ds['lat'].values
    lons = Tmin_ds['lon'].values
    lat_indices = np.where((lats >= min_lat) & (lats <= max_lat))[0]
    lon_indices = np.where((lons >= min_lon) & (lons <= max_lon))[0]
    Tmin = Tmin_ds['Temperature_Air_2m_Min_24h'].values[0, lat_indices.min():lat_indices.max()+1, lon_indices.min():lon_indices.max()+1]
    Tmax = Tmax_ds['Temperature_Air_2m_Max_24h'].values[0, lat_indices.min():lat_indices.max()+1, lon_indices.min():lon_indices.max()+1]
    return Tmin, Tmax, lats, lons, lat_indices, lon_indices

def read_sr(year, month, day, sr_dir, min_lat, max_lat, min_lon, max_lon):
    # Variables in SR file: ['Solar_Radiation_Flux', 'time', 'lon', 'lat']
    target_date = f"{year}{month.zfill(2)}{str(day).zfill(2)}"
    sr_file = os.path.join(sr_dir, f"Solar-Radiation-Flux_C3S-glob-agric_AgERA5_{target_date}_final-v1.1.area-subset.-20.15.-35.35.nc")
    if not os.path.exists(sr_file):
        print(f"Missing SR file for {year}-{month.zfill(2)}-{str(day).zfill(2)}. Skipping this day.")
        return None
    sr_ds = xr.open_dataset(sr_file)
    lats = sr_ds['lat'].values
    lons = sr_ds['lon'].values
    lat_indices = np.where((lats >= min_lat) & (lats <= max_lat))[0]
    lon_indices = np.where((lons >= min_lon) & (lons <= max_lon))[0]
    SR = sr_ds['Solar_Radiation_Flux'].values[0, lat_indices.min():lat_indices.max()+1, lon_indices.min():lon_indices.max()+1]
    return SR, lats, lons, lat_indices, lon_indices

def read_ws(year, month, day, ws_dir, min_lat, max_lat, min_lon, max_lon):
    # Variables in WS file: ['Wind_Speed_10m_Mean', 'time', 'lon', 'lat']
    target_date = f"{year}{month.zfill(2)}{str(day).zfill(2)}"
    ws_file = os.path.join(ws_dir, f"Wind-Speed-10m-Mean_C3S-glob-agric_AgERA5_{target_date}_final-v1.1.area-subset.-20.15.-35.35.nc")
    if not os.path.exists(ws_file):
        print(f"Missing WS file for {year}-{month.zfill(2)}-{str(day).zfill(2)}. Skipping this day.")
        return None
    ws_ds = xr.open_dataset(ws_file)
    lats = ws_ds['lat'].values
    lons = ws_ds['lon'].values
    lat_indices = np.where((lats >= min_lat) & (lats <= max_lat))[0]
    lon_indices = np.where((lons >= min_lon) & (lons <= max_lon))[0]
    WS = ws_ds['Wind_Speed_10m_Mean'].values[0, lat_indices.min():lat_indices.max()+1, lon_indices.min():lon_indices.max()+1]
    return WS, lats, lons, lat_indices, lon_indices

def read_pr(year, month, day, pr_dir, min_lat, max_lat, min_lon, max_lon):
    # Variables in Pr file: ['Precipitation_Flux', 'time', 'lon', 'lat']
    target_date = f"{year}{month.zfill(2)}{str(day).zfill(2)}"
    pr_file = os.path.join(pr_dir, f"Precipitation-Flux_C3S-glob-agric_AgERA5_{target_date}_final-v1.1.area-subset.-20.15.-35.35.nc")
    if not os.path.exists(pr_file):
        print(f"Missing Pr file for {year}-{month.zfill(2)}-{str(day).zfill(2)}. Skipping this day.")
        return None
    pr_ds = xr.open_dataset(pr_file)
    lats = pr_ds['lat'].values
    lons = pr_ds['lon'].values
    lat_indices = np.where((lats >= min_lat) & (lats <= max_lat))[0]
    lon_indices = np.where((lons >= min_lon) & (lons <= max_lon))[0]
    Pr = pr_ds['Precipitation_Flux'].values[0, lat_indices.min():lat_indices.max()+1, lon_indices.min():lon_indices.max()+1]
    return Pr, lats, lons, lat_indices, lon_indices

def calculate_eto(z, lats, lons, Tmin, Tmax, VP, SR, WS):
    # Convert lats and lons to xarray DataArrays if they are not already
    if not isinstance(lats, xr.DataArray):
        lats = xr.DataArray(lats, dims="lat", coords={"lat": lats})
    if not isinstance(lons, xr.DataArray):
        lons = xr.DataArray(lons, dims="lon", coords={"lon": lons})

    # Ensure lats and lons are broadcast to the same shape as the grid
    lat_grid, lon_grid = xr.broadcast(lats, lons)
    # Define constants
    G = 0                 # Soil heat flux density (MJ/m²/day)
    Gsc = 0.0820          # Solar constant (MJ/m²/min)
    sigma = 4.903e-9      # Stefan-Boltzmann constant (MJ/m²/K⁴/day)
    cp = 1.013e-3         # Specific heat capacity of air (MJ/kg/°C)
    epsilon = 0.622       # Ratio of molecular weight of water vapour/dry air
    Lv = 2.45             # Latent heat of vaporization (MJ/kg)
    alpha = 0.23          # Albedo for hypothetical grass reference crop       
    # Unit conversions
    phi = np.radians(lat_grid)  # Latitude from ° to rad
    TminC = Tmin - 273.15   # Convert from K to °C
    TmaxC = Tmax - 273.15   # Convert from K to °C
    Rs = SR * 1e-6          # Convert from J/m²/day to MJ/m²/day
    U2 = WS * 0.748         # Convert wind speed from 10 m to 2 m
    ea = VP * 0.1           # Convert from hPa to kPa
    # Mean temperature
    Tmean = (Tmin + Tmax) / 2      # (K)
    TmeanC = (TminC + TmaxC) / 2   # (°C) 
    # Atmospheric pressure (kPa)
    P = 101.3 * ((Tmean - 0.0065 * z) / Tmean) ** 5.26   # use Tmean in K
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
    Rso = (0.75 + 2e-5 * z) * Ra
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
    return ETo

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

# Retreive elevation data
z_file = os.path.join(dem_dir, "AgERA5_GridPoint_Elevations.csv")
elevation_data = pd.read_csv(z_file)
z = None  # Initialize elevation grid

for day in range(1, days_in_month + 1):
    
    # Calculate RHmin and RHmax and save these to individual files
    # RHmin and RHmax are not used in ETo calculation since VP is available
    # However, these values are required for the database
    RHmin, RHmax, lats, lons, lat_indices, lon_indices, ds = calculate_rhmin_rhmax(
        year, month, day, rh_dir, min_lat, max_lat, min_lon, max_lon)
    if RHmin is None or RHmax is None:
        continue

    output_ds_min = xr.Dataset(
        {"RHmin": (["lat", "lon"], RHmin)},
        coords={"lat": lats[lat_indices.min():lat_indices.max()+1], "lon": lons[lon_indices.min():lon_indices.max()+1]},
        attrs=ds.attrs
    )

    output_ds_max = xr.Dataset(
        {"RHmax": (["lat", "lon"], RHmax)},
        coords={"lat": lats[lat_indices.min():lat_indices.max()+1], "lon": lons[lon_indices.min():lon_indices.max()+1]},
        attrs=ds.attrs
    )

    target_date = f"{year}{month.zfill(2)}{str(day).zfill(2)}"
    rhmin_filename = os.path.join(rh_dir, f"RHmin-2m_C3S-glob-agric_AgERA5_{target_date}_final-v1.1.area-subset.-20.15.-35.35.nc")
    rhmax_filename = os.path.join(rh_dir, f"RHmax-2m_C3S-glob-agric_AgERA5_{target_date}_final-v1.1.area-subset.-20.15.-35.35.nc")

    output_ds_min.to_netcdf(rhmin_filename, engine="netcdf4")
    output_ds_max.to_netcdf(rhmax_filename, engine="netcdf4")

    # Read Vapour Pressure (VP)
    VP, lats_vp, lons_vp, lat_indices_vp, lon_indices_vp = read_vp(
        year, month, day, vp_dir, min_lat, max_lat, min_lon, max_lon)
    if VP is None:
        continue

    # Read Tmin and Tmax and convert to degrees Celsius
    Tmin, Tmax, lats_temp, lons_temp, lat_indices_temp, lon_indices_temp = read_tmin_tmax(
    year, month, day, temp_dir, min_lat, max_lat, min_lon, max_lon)
    if Tmin is None or Tmax is None:
        continue
    
    # Read Solar Radiation (SR)
    SR, lats_sr, lons_sr, lat_indices_sr, lon_indices_sr = read_sr(
        year, month, day, sr_dir, min_lat, max_lat, min_lon, max_lon)
    if SR is None:
        continue

    # Read Wind Speed (WS)
    WS, lats_ws, lons_ws, lat_indices_ws, lon_indices_ws = read_ws(
        year, month, day, ws_dir, min_lat, max_lat, min_lon, max_lon)
    if WS is None:
        continue

    # Read Precipitation (Pr)
    Pr, lats_pr, lons_pr, lat_indices_pr, lon_indices_pr = read_pr(
        year, month, day, pr_dir, min_lat, max_lat, min_lon, max_lon)
    if Pr is None:
        continue
    
    # Extract elevation data for each gridpoint (once only)
    if z is None:  # Check if z has already been calculated
        z = np.full((len(lats), len(lons)), np.nan)  # Initialize z with NaN
        z = get_elevation(z, lats, lons, elevation_data)  # Calculate z for the first time
    
    # Calculate ETo
    #ETo = calculate_eto(z, lats, lons, Tmin, Tmax, VP, SR, WS)
    ETo = calculate_eto(z, lats, lons, Tmin, Tmax, VP, SR, WS)
    # Create dataset and save to NetCDF
    eto_dataset = xr.Dataset(
    {"ETo": (["lat", "lon"], ETo.data if isinstance(ETo, xr.DataArray) else ETo)},
    coords={
        "lat": lats.data if isinstance(lats, xr.DataArray) else lats,
        "lon": lons.data if isinstance(lons, xr.DataArray) else lons
            },
    attrs={"description": "Reference Evapotranspiration (PM-FAO56)"}
    )
    eto_filename = os.path.join(eto_dir, f"ETo_C3S-glob-agric_AgERA5_{year}{month.zfill(2)}{str(day).zfill(2)}_final-v1.1.area-subset.-20.15.-35.35.nc")
    eto_dataset.to_netcdf(eto_filename, engine="netcdf4")
    
    # Message that shows progress
    print(f"Processed {year}-{month.zfill(2)}-{str(day).zfill(2)} successfully.")
