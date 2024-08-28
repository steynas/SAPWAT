# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 16:00:38 2024

@author: SteynAS
"""

import cdsapi
import os
import calendar

def download_RH(year, month, base_dir):
    # Define the download directory based on the base_dir
    download_directory = os.path.join(base_dir, "RH")
    if not os.path.exists(download_directory):
        os.makedirs(download_directory)

    # Determine the number of days in the selected month
    days_in_month = calendar.monthrange(int(year), int(month))[1]

    # Create the list of days based on the month
    days = [str(day).zfill(2) for day in range(1, days_in_month + 1)]

    # Define the CDS API request
    dataset = "sis-agrometeorological-indicators"
    request = {
        'variable': '2m_relative_humidity',
        'year': [year],
        'month': [month],
        'day': days,
        'time': ['06_00', '09_00', '12_00', '15_00', '18_00'],
        'version': '1_1',
        'area': [-20, 15, -35, 35]
    }

    # Create a more descriptive file name based on the year and month
    month_name = calendar.month_name[int(month)]
    output_file = os.path.join(download_directory, f"AgERA5_RH_{month_name}_{year}.zip")

    # Download the data using the CDS API
    client = cdsapi.Client()
    client.retrieve(dataset, request).download(output_file)

    print(f"RH data successfully downloaded to {output_file}")
