#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=6400
#SBATCH --time=01:00:00
#SBATCH --job-name=my_python_job
#SBATCH --output=my_python_job_%j.out

# Load the required modules
module load Python

# Set up a temporary directory for the virtual environment
temp_dir="/home"
mkdir -p $temp_dir

# Create and activate the virtual environment
python -m venv $temp_dir/venv
source $temp_dir/venv/bin/activate

# Install required packages in the virtual environment
pip install pandas

# Change to the directory where your Python script and data are located
cd /home

# Run the Python script using the appropriate command
python abhishek_python_script.py

# Deactivate and remove the virtual environment
deactivate
rm -rf $temp_dir/venv

