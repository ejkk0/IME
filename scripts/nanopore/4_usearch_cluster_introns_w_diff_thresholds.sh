#!/bin/bash

#SBATCH --partition sched_mit_cburge
#SBATCH --job-name=usearch_cluster_introns
#SBATCH --output=/home/ejkk/Cluster/IME/outfiles/usearch_cluster_filtered_introns_diff_thresholds.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=128000

echo "hola"

input_file="/pool001/ejkk/RIL2/filtered_introns.fa"
output_dir="/pool001/ejkk/RIL2/"

# Array of identity thresholds
identity_thresholds=(0.5 0.6 0.7 0.8 0.9)

for threshold in "${identity_thresholds[@]}"; do
    # Define output filenames based on the threshold
    consout_file="${output_dir}f_cons_${threshold//./}.fa"
    msaout_file="${output_dir}MSAs/f_msa_${threshold//./}.fa"
    uc_file="${output_dir}f_clusters_${threshold//./}.uc"

    # Run usearch command with the current threshold
    usearch -cluster_fast "$input_file" -id "$threshold" -consout "$consout_file" -msaout "$msaout_file" -uc "$uc_file"

    echo "Clustering for threshold $threshold done!"

    # Merge MSAs into one file
    tail -n +1 "${output_dir}MSAs/f_msa_${threshold//./}.fa"* > "${output_dir}MSAs/f_msa_${threshold//./}_all.fa"

    echo "Merging MSAs for threshold $threshold done!"

    # Delete individual MSA files
    rm "${output_dir}MSAs/f_msa_${threshold//./}.fa"*
    echo "Deleting individual MSA files for threshold $threshold done!"
done

echo "all done :)"
