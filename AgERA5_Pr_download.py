# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 22:03:26 2024

@author: SteynAS
"""

import cdsapi
import os
import calendar

def download_Pr(year, month, base_dir):
    # Define the download directory based on the base_dir
    download_directory = os.path.join(base_dir, "Pr")
    if not os.path.exists(download_directory):
        os.makedirs(download_directory)

    # Determine the number of days in the selected month
    days_in_month = calendar.monthrange(int(year), int(month))[1]

    # Create the list of days based on the month
    days = [str(day).zfill(2) for day in range(1, days_in_month + 1)]

    # Define the CDS API request for daily total precipitation flux
    dataset = "sis-agrometeorological-indicators"
    request = {
        'variable': 'precipitation_flux',
        'year': [year],
        'month': [month],
        'day': days,
        'version': '1_1',
        'area': [-20, 15, -35, 35]
    }

    # Create a more descriptive file name based on the year and month
    month_name = calendar.month_name[int(month)]
    output_file = os.path.join(download_directory, f"AgERA5_Pr_{month_name}_{year}.zip")

    # Download the data using the CDS API
    client = cdsapi.Client()
    client.retrieve(dataset, request).download(output_file)

    print(f"Precipitation Flux data successfully downloaded to {output_file}")
