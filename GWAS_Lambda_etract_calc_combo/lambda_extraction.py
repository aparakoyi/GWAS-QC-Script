#!/usr/bin/env Python3.10
# lambda_extraction.py

import os
import gwaslab
import numpy as np
import scipy as sp


# Function to extract lambda value for a GWAS file
def extract_lambda(f_path):
    try:
        # Load the GWAS summary statistics file
        df = gwaslab.Sumstats(f_path, snpid="SNP_ID", chrom="chrom", pos="pos",
                              ea='ea', ref="ref", alt="alt", eaf="af", p="pval", OR="or")

        # calculate z score
        z_score = sp.stats.norm.ppf(df.data["P"] / 2)

        # Calculate lambda
        lb = round(np.median(z_score ** 2) / 0.454, 3)

        return lb

    except Exception as e:
        print(f"Error processing {f_path}: {e}")
        return None


# Directory containing GWAS summary statistics files
input_directory = "/Volumes/AbigailHD/MVP_R4.1000G_AGR.GIA.PheCodes_CirculatorySystem_batch1/EUR"
output_file = "/Volumes/AbigailHD/lambda_Val_files/CirculatorySystem_batch1_EUR_lambda_values_table.txt"

# Open output file for writing
with open(output_file, "w") as out:
    # Write header
    out.write("Filename\tLambda\n")

    # Process each file in the directory
    for file_name in os.listdir(input_directory):
        file_path = os.path.join(input_directory, file_name)

        # Skip if not a file
        if not os.path.isfile(file_path):
            continue

        # Extract lambda value
        lambda_value = extract_lambda(file_path)

        # Write result to the output file
        if lambda_value is not None:
            out.write(f"{file_name}\t{lambda_value:.4f}\n")

print(f"Lambda values have been saved to {output_file}.")
