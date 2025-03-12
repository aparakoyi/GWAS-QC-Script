# GWAS-QC

## About The Project
This project is intended to be used for the purpose of quality checking and standardizing GWAS files.
Most of the functionalities within this project are built around the [gwaslab](https://github.com/Cloufield/gwaslab/) package. 

## Getting Started
In preparation to run this script ensure that all gwas files are placed in one data folder and there is an empty output folder 
present within said data folder. Following the creation of the data folder complete the yaml file for the respective data.
Record the paths to the gwas file directroy, the outfile directory, and the yaml config file as these will be used to run the 
script. 
### Prerequisites 
Below are the necessary software packages, versions, and files necessary to run the script. 
* Gwaslab
````commandline
pip install gwaslab
````
* Argparse
```commandline
pip install argparse
```
* Yaml
```commandline
pip install yaml
```
* Python == 3.10, this script was made using this specific python package so creation of a conda environment is recommanded. 
* Config yaml file 

If running this script on an HPC i.e. Eristwo, below are the additional packages to run on Eristwo to run this script.
* Miniforge3
```commandline
module load minifogre3
```
If using a conda environment 
* Conda environment
```commandline
conda activate <environment_name>
```
* Reference allele and VCF folder 
## Usage
The QC script takes in three command line inputs
Within the folder that holds both the GWAS_QC script and GWAS_QC script yaml config file.
* -i, --infile : path to GWAS file(s)
* -c, --config : path to yaml config file 
* -o, --outfile : path and name to quality checked file(s)

Below is an example of the command line code one would run to run the script. 
```commandline
New_GWAS_QC_Code.py -i path_of_GWAS_files_of_interest -c path_of_GWAS_QC_config.yaml -o path_of_GWAS_ouput_folder
```
If running this script using an HPC cluster you would run it using a bash script. Below is an example bash script for 
submitting an array task containing files to be QC'd.
```bash
#!/bin/bash

#SBATCH --job-name= your_job_name 

#SBATCH --partition=normal

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=40

#SBATCH --array= # - ## (your_range_of_array_tasks

#SBATCH --error=log%j.err

# specify array config file
a_config= path_to_array_job_list_file.txt

# extract the file name for the current $SLURM_ARRAY_TASK_ID
file_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $a_config)

# print to stdout the current slurm task log
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the file name is ${file_name}." >> /dev/stdout

# Run python program
python New_GWAS_QC_Code.py -i path_of_GWAS_files_of_interest/${file_name} 
-c path_of_GWAS_QC_config.yaml -o path_of_GWAS_ouput_folder/${file_name}
```
## TroubleShooting 
* Ensure that you installed the latest version of gwaslab and that you are using a python version less than or equal to 
3.10. 

This QC script is divided into 4 main sections (Basic check, liftover, stats check, and Harmonization). Below are a few
general troubleshooting issues that could occur in these 4 sections of the code. 
* Basic Check 
  * The GWAS lab package requires at least a SNIPID or rsID column to be included in any dataset. If either of these are
    not included or set when loading in the data the file will not be able to be loaded in. This is the importance of the 
    yaml file as it allows for the specification of the necessary columns that the GWASlab package requires. 
  * The basic check function does not change or remove anything but rather documents everything check in the log file. 
    Ensure you always check the log file for checks its made. 
* Lift over 
  * If the script you are running is already in build 38 this script should not do lift over but in teh event that liftover
  occurs on a file that is in build 38 already the CHR and POS columns will be completely wrong and change to build 99. 
  * This liftover step requires API calls so ensure whatever device you are suing is connected to the internet for this 
  step to run. 
* Stats Check 
  * In order for any of the stats check to run we require specific columns for each. 
    * For mlog10P --> P, we require the dataset to have a mlog10p column
    * For OR --> Beta, we require the dataset to have a OR column 
    * For maf calculation --> we require the dataset to have an EAF column 
  * We cannot calculate the EAF or MAF without at least one of these columns present using this script. 
* Harmonization
  * Ensure that the directory you are working in holds the Ref allele and vcf file folder. 
  * This step takes is the longest step in the script. 
  * Checking of palindromic SNPs and indels occurs here, but we do not make any changes to these. 
  

## Acknowledgements
 * Abigail Parakoyi
 * United States Veteran Assistance Department Drug Discovery Group 
 * GWASLab preprint: He, Y., Koido, M., Shimmori, Y., Kamatani, Y. (2023). GWASLab: a Python package for processing and 
    visualizing GWAS summary statistics. Preprint at Jxiv, 2023-5. https://doi.org/10.51094/jxiv.370