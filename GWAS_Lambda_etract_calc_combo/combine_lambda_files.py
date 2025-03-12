#!/usr/bin/env Python3.10
# combine_lambda_files.py

import os
import pandas as pd


# combine data frames
def combine_df(input_dir, output_file):
    combined_df = []

    # Process each file in the directory
    for file_name in os.listdir(input_dir):
        file_path = os.path.join(input_dir, file_name)
    # Skip if not a file
        if not os.path.isfile(file_path):
            continue
        try:
            # Read the file into a DataFrame
            df = pd.read_csv(file_path, sep='\t')

            # Append to combined data
            combined_df.append(df)
        except Exception as e:
            print(f"Error processing {file_path}: {e}")

    # Concatenate all dataframes into one
    if combined_df:
        combined_df = pd.concat(combined_df, ignore_index=True)
        # sort data by lambda
        if 'Lambda' in combined_df.columns:
            combined_df = combined_df.sort_values(by='Lambda', ascending=False)
        # Save the combined table
        combined_df.to_csv(output_file, sep='\t', index=False)
        print(f"Combined data table saved to {output_file}")
    else:
        print("No data tables found to combine.")


input_dir = "/Volumes/AbigailHD/lambda_Val_files"
output_file = "Combined_lambda_table.txt"

# combine tables
combine_df(input_dir, output_file)
