# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 12:39:19 2024

@author: SteynAS
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 21:22:07 2024

@author: SteynAS
"""

from AgERA5_T_download import download_T
from AgERA5_VP_download import download_VP
from AgERA5_RH_download import download_RH
from AgERA5_WS_download import download_WS
from AgERA5_SR_download import download_SR
from AgERA5_Pr_download import download_Pr

def main():
    # Get user input for the start year and end year
    start_year = int(input("Enter the start year (e.g., 1979): "))
    end_year = int(input("Enter the end year (e.g., 2024): "))

    # Base directory
    base_dir = "/home/steynas/AgERA5/home/steynas/AgERA5"
    
    # Loop through each year and each month
    for year in range(start_year, end_year + 1):
        for month in range(1, 13):
            month_str = str(month).zfill(2)  # Format month as 01, 02, ..., 12
            # Show progress
            print("Processing: ",year , month_str)
            
            # Download Tmin and Tmax data
            download_T(str(year), month_str, base_dir)
    
            # Download Vapour Pressure data
            download_VP(str(year), month_str, base_dir)
    
            # Download Relative Humidity data
            download_RH(str(year), month_str, base_dir)
    
            # Download Wind Speed data
            download_WS(str(year), month_str, base_dir)
    
            # Download Solar Radiation Flux data
            download_SR(str(year), month_str, base_dir)
    
            # Download Precipitation Flux data
            download_Pr(str(year), month_str, base_dir)

if __name__ == "__main__":
    main()
