# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 09:00:20 2024

This scripts reads the daily station data from Excel (combined) = stn
and the climate database reanalysis data from CSV (nearest grid point) =cdb
and merge them into a single Excel file for statistical analysis
Just change the name of the input and output file & path, e.g.:
input1 = Glen_Combined_1984_2023.xlsx
input2 = ClimateTimeSeries_AgERA5_-28.9_26.3.csv
output = Glen_Merged_stn_cdb.xlsx

@author: SteynAS
"""

import os
import pandas as pd

# Define the directory and file names
stn_dir = r'C:\StnData'
cdb_dir = r'C:\AgERA5\ClimDataBase'
stn_file = 'Glen_Combined_1984_2023.xlsx'
cdb_file = 'ClimateTimeSeries_AgERA5_-28.9_26.3.csv'
out_file = 'Glen_Merged_stn_cdb.xlsx'
output_file_path = os.path.join(stn_dir, out_file)

# Paths to the files
file_path_cdb = os.path.join(cdb_dir, cdb_file)
file_path_stn = os.path.join(stn_dir, stn_file)

# Load the CSV file into a DataFrame
df_cdb = pd.read_csv(file_path_cdb, parse_dates=['Date'])

# Load the Excel file into a DataFrame
df_stn = pd.read_excel(file_path_stn, parse_dates=['Date'])

# Display basic information about the CSV dataset
print("\nBasic Information about the Climate Time Series CSV dataset:")
df_cdb.info()

# Display the first few rows of the CSV dataset
print("\nFirst 5 rows of the Climate Time Series CSV dataset:")
print(df_cdb.head())

# Display basic information about the Excel dataset
print("\nBasic Information about the Glen Combined Excel dataset:")
df_stn.info()

# Display the first few rows of the Excel dataset
print("\nFirst 5 rows of the Glen Combined Excel dataset:")
print(df_stn.head())

# Define the complete date range for alignment
date_range = pd.date_range(start='1979-01-01', end='2023-12-31')

# Set 'Date' as the index for both dataframes
df_cdb.set_index('Date', inplace=True)
df_stn.set_index('Date', inplace=True)

# Reindex the cdb dataframe to ensure it covers the complete date range
# stn dataframe is kept as-is since we only want rows where data exists.
df_cdb = df_cdb.reindex(date_range)

# Reset the index to have 'Date' as a column again
df_cdb.reset_index(inplace=True)
df_stn.reset_index(inplace=True)

# Rename the index column back to 'Date'
df_cdb.rename(columns={'index': 'Date'}, inplace=True)
df_stn.rename(columns={'index': 'Date'}, inplace=True)

# Rename columns as specified
# Station (stn) data renaming
stn_rename_mapping = {
    'Tx': 'stn_Tmax',
    'Tn': 'stn_Tmin',
    'RHx': 'stn_RHmax',
    'RHn': 'stn_RHmin',
    'U2': 'stn_WSmean',
    'Rs est.': 'stn_SR',
    'Rain': 'stn_Pr',
    'PM ET0': 'stn_ETo'
}
df_stn.rename(columns=stn_rename_mapping, inplace=True)

# Climate database (cdb) data renaming
cdb_rename_mapping = {
    'Tmax': 'cdb_Tmax',
    'Tmin': 'cdb_Tmin',
    'RHmax': 'cdb_RHmax',
    'RHmin': 'cdb_RHmin',
    'WSmean': 'cdb_WSmean',
    'SR': 'cdb_SR',
    'Pr': 'cdb_Pr',
    'ETo': 'cdb_ETo'
}
df_cdb.rename(columns=cdb_rename_mapping, inplace=True)

# Merge the dataframes on 'Date' for comparison
df_merged = pd.merge(df_stn, df_cdb, on='Date', how='left', suffixes=('_stn', '_cdb'))

# Filter the merged dataframe to keep only the rows where stn data exists (non-null 'stn_Tmax')
df_merged = df_merged[df_merged['stn_Tmax'].notnull()]

# Select only the columns to include in the output
columns_to_include = ['Date'] + list(stn_rename_mapping.values()) + list(cdb_rename_mapping.values())
df_merged = df_merged[columns_to_include]

# Save the filtered merged dataset to the specified output file
df_merged.to_excel(output_file_path, index=False)
print(f"\nMerged comparison data saved to: {output_file_path}")

# Display the first few rows of the final output for verification
print("\nFirst 5 rows of the final output:")
print(df_merged.head())
