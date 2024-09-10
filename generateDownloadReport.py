# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 13:39:43 2024

@author: SteynAS
"""
import os
import csv
from datetime import datetime

# Function to convert month number to full month name
def month_name(month):
    months = ["January", "February", "March", "April", "May", "June", 
              "July", "August", "September", "October", "November", "December"]
    return months[month - 1]

# Function to check if all required files are present for a given year and month
def check_files(year, month, directories):
    month_str = month_name(month)  # Convert month number to month name
    
    # Define the expected file names for each directory
    expected_files = [
        f"AgERA5_Temps_{month_str}_{year}.zip",        # Temps (Tmin, Tmax)
        f"AgERA5_VP_{month_str}_{year}.zip",           # Vapour Pressure
        f"AgERA5_RH_{month_str}_{year}.zip",           # Relative Humidity
        f"AgERA5_WS_{month_str}_{year}.zip",           # Wind Speed
        f"AgERA5_SR_{month_str}_{year}.zip",           # Solar Radiation
        f"AgERA5_Pr_{month_str}_{year}.zip"            # Precipitation
    ]
    
    # Check each directory for the corresponding expected file
    all_files_present = True
    for dir_path, expected_file in zip(directories, expected_files):
        file_path = os.path.join(dir_path, expected_file)
        
        # Print the file path for debugging
        print(f"Checking: {file_path}")
        
        if not os.path.exists(file_path):
            print(f"Missing: {file_path}")  # Print missing file for debugging
            all_files_present = False
    
    return all_files_present  # Return True if all required files are present

def generate_report(start_year, end_year, directories, base_dir):
    # Initialize the headers and data rows
    headers = ["Year", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    rows = []

    # Loop through each year and month to check for file presence
    for year in range(start_year, end_year + 1):
        row = [year]  # Start the row with the year
        for month in range(1, 12 + 1):  # Loop through each month
            if check_files(year, month, directories):
                row.append("Yes")
            else:
                row.append("")  # Leave blank if any file is missing
        rows.append(row)

    # Generate the file name with the current date
    current_date = datetime.now().strftime("%Y%m%d")
    output_csv = os.path.join(base_dir, f"DownloadReport_{current_date}.csv")
    
    # Write the data to a CSV file
    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(headers)  # Write the headers
        writer.writerows(rows)  # Write the data rows

    print(f"Report generated: {output_csv}")

# Example usage
base_dir = "/home/steynas/AgERA5/"
directories = [
    os.path.join(base_dir, "Temps"),  # Directory for Temperature (Tmin, Tmax)
    os.path.join(base_dir, "VP"),     # Directory for Vapour Pressure
    os.path.join(base_dir, "RH"),     # Directory for Relative Humidity (06h, 09h, 12h, 15h, 18h)
    os.path.join(base_dir, "WS"),     # Directory for Wind Speed
    os.path.join(base_dir, "SR"),     # Directory for Solar Radiation
    os.path.join(base_dir, "Pr")      # Directory for Precipitation
]

# Run the report generation and save it to the base directory on the UFS cluster
generate_report(1979, 2024, directories, base_dir)
