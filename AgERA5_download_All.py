# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 21:22:07 2024

@author: SteynAS
"""

from AgERA5_T_download import download_T
from AgERA5_RH_download import download_RH
from AgERA5_WS_download import download_WS
from AgERA5_SR_download import download_SR
from AgERA5_Pr_download import download_Pr

def main():
    # Get user input for the year and month
    year = input("Enter the year (e.g., 1979): ")
    month = input("Enter the month (01-12): ")
    
    # Base directory
    base_dir = "C:/AgERA5"
    
    # Download Tmin and Tmax data
    download_T(year, month, base_dir)
    
    # Download Relative Humidity data
    download_RH(year, month, base_dir)
    
    # Download Wind Speed data
    download_WS(year, month, base_dir)
    
    # Download Solar Radiation Flux data
    download_SR(year, month, base_dir)
    
    # Download Precipitation Flux data
    download_Pr(year, month, base_dir)

if __name__ == "__main__":
    main()