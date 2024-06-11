# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 12:31:05 2024

@author: steynas
"""

import pandas as pd
import os

# Define input and output directories
input_dir = r"C:\Data\ETo stations\Before Python"
output_dir = r"C:\Data\ETo stations\Processed"

# Define input and output file names
input_file = os.path.join(input_dir, "Tmax_Cradock.xlsx")
output_file = os.path.join(output_dir, "Tmax_Cradock_processed.xlsx")

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

# Debug: Print the unique years to check if 1985 is included
print("Unique years in combined data before melting:", df_combined['Year'].unique())

# Melt the DataFrame to have a 'Month' column and a 'Tmax' column
df_melted = df_combined.melt(id_vars=['Station', 'Year', 'Day'],
                             value_vars=['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'],
                             var_name='Month',
                             value_name='Tmax')

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

# Debug: Print unique years to check if 1985 is included
print("Unique years in expanded data:", df_expanded['Year'].unique())

# Fill the 'Station' column with the correct station name
df_expanded['Station'] = df_expanded['Station'].ffill()

# Write the processed data to a new Excel file with a single header
with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
    df_expanded[['Station', 'Year', 'Month', 'Day', 'Tmax']].to_excel(writer, index=False, sheet_name='Tmax')

print(f"Processed data has been written to {output_file}")