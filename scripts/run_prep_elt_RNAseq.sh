#!/bin/bash

#SBATCH --partition sched_mit_cburge
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=512000
#SBATCH --mail-type=END
#SBATCH --mail-user=ejkk@mit.edu

# Job name and output file name based on command-line args
#SBATCH --job-name=preprocess_RNAseq_nbc_3way_spl_filter
#SBATCH --output=/home/ejkk/Cluster/IME/outfiles/%x_%j.out

# Check if the number of provided arguments is less than 2 and set defaults
if [ $# -lt 2 ]; then
  echo "Not all arguments provided, using default values for missing arguments."
  # Default values if not all arguments are provided
  ARG1=${1:-'test'} # Default for var1; can also be 'RIL2' or 'TRIS'
  ARG2=${2:-1}   # Default for var2; can also be 2
else
  ARG1=$1
  ARG2=$2
fi

# Run the Python script
python preprocess_element_RNAseq_no_barcode_filter.py $ARG1 $ARG2
