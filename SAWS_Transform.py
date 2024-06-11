# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 15:56:35 2024

@author: steynas
"""

import pandas as pd
import os
import sys

# Define input and output directories
input_dir = r"C:\Data\ETo stations\Before Python"
output_dir = r"C:\Data\ETo stations\Processed"

# Prompt user to choose the station and climatic element
station = input("Enter the station name (e.g., Cradock): ")

print("Choose a climatic element from the following options:")
climatic_elements = ["Tmin", "Tmax", "Rain", "RH_08h00", "RH_14h00", "RH_20h00", 
                     "Wind_08h00", "Wind_14h00", "Wind_20h00", "Rad"]

for i, element in enumerate(climatic_elements, 1):
    print(f"{i}. {element}")

element_choice = int(input("Enter the number corresponding to the climatic element: "))
value_name = climatic_elements[element_choice - 1]

# Define input and output file names based on user input
input_file = os.path.join(input_dir, f"{value_name}_{station}.xlsx")
output_file = os.path.join(output_dir, f"{value_name}_{station}_processed.xlsx")

# Check if the input file exists
if not os.path.exists(input_file):
    print(f"Error: The input file '{input_file}' does not exist. Please check the station name and climatic element, and ensure the file is in the specified path.")
    sys.exit(1)

# Read the Excel file
df = pd.read_excel(input_file)

# Strip any leading/trailing whitespace from the column names
df.columns = df.columns.str.strip()

# Identify header rows
header_rows = df[df['Day'] == 'Day'].index.tolist()

# Remove rows that are headers repeating within the data
df = df[~df.index.isin(header_rows)]

# Process the data
processed_data = []
start_idx = 0
for i in range(len(header_rows)):
    end_idx = header_rows[i] if i < len(header_rows) else len(df)
    year_data = df.iloc[start_idx:end_idx].copy()
    year_data.columns = ['Station', 'Year', 'Day', 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
    processed_data.append(year_data)
    start_idx = end_idx + 1

# Combine all processed data
df_combined = pd.concat(processed_data, ignore_index=True)

# Replace -9999 with NaN for missing values
df_combined.replace(-9999, pd.NA, inplace=True)

# Convert columns to appropriate types
df_combined['Year'] = df_combined['Year'].astype(int)
df_combined['Day'] = df_combined['Day'].astype(int)

# Print the selected station and climatic element
min_year = df_combined['Year'].min()
max_year = df_combined['Year'].max()
print(f"Station chosen: {station}")
print(f"Climatic element selected: {value_name}")
print(f"Start date: {min_year}")
print(f"End date: {max_year}")

# Melt the DataFrame to have a 'Month' column and a 'value_name' column
df_melted = df_combined.melt(id_vars=['Station', 'Year', 'Day'],
                             value_vars=['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'],
                             var_name='Month',
                             value_name=value_name)

# Map month names to numbers
month_mapping = {'JAN': 1, 'FEB': 2, 'MAR': 3, 'APR': 4, 'MAY': 5, 'JUN': 6,
                 'JUL': 7, 'AUG': 8, 'SEP': 9, 'OCT': 10, 'NOV': 11, 'DEC': 12}
df_melted['Month'] = df_melted['Month'].map(month_mapping)

# Ensure the 'Month' column is of the same type in both dataframes
df_melted['Month'] = df_melted['Month'].astype(int)

# Generate a complete date range from the minimum to maximum year in the original data
min_year = df_combined['Year'].min()
max_year = df_combined['Year'].max()

# Generate all valid dates between min_year and max_year
all_dates = pd.date_range(start=f'{min_year}-01-01', end=f'{max_year}-12-31', freq='D')
all_dates_df = pd.DataFrame({
    'Year': all_dates.year,
    'Month': all_dates.month,
    'Day': all_dates.day
})

# Merge the original data with the complete date range
df_expanded = pd.merge(all_dates_df, df_melted, how='left', on=['Year', 'Month', 'Day'])

# Fill the 'Station' column with the correct station name
df_expanded['Station'] = df_expanded['Station'].ffill()

# Replace missing values with -9999
missing_values_count = df_expanded.isna().sum().sum()
df_expanded.fillna(-9999, inplace=True)

# Print the number of missing values
print(f"Missing values: {missing_values_count}")
print("Missing values indicated by -9999")

# Write the processed data to a new Excel file with a single header
with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
    df_expanded[['Station', 'Year', 'Month', 'Day', value_name]].to_excel(writer, index=False, sheet_name=value_name)

print(f"Processed data has been written to {output_file}")
