# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 15:38:28 2025

@author: steynas
"""

import os
import pandas as pd

def repair_csv_dates(directory="/home/steynas/AgERA5/ClimDataBase"):
    """
    Fixes all CSV files in the given directory by ensuring the Date column is in YYYYMMDD format.
    Only modifies files where the Date format was wrongly converted to YYYY-MM-DD.
    
    :param directory: Path to the directory containing CSV files.
    """
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            file_path = os.path.join(directory, filename)
            
            try:
                df = pd.read_csv(file_path, dtype={"Date": str})
                
                # Check if the Date format is incorrect (contains dashes)
                if df["Date"].str.contains("-").any():
                    print(f"Repairing date format in: {filename}")
                    
                    # Convert Date column back to YYYYMMDD format
                    df["Date"] = pd.to_datetime(df["Date"], errors='coerce').dt.strftime('%Y%m%d')
                    
                    # Save back to the same file
                    df.to_csv(file_path, index=False)
                    print(f"Repaired: {filename}")
            except Exception as e:
                print(f"Error processing {filename}: {e}")

if __name__ == "__main__":
    repair_csv_dates()
