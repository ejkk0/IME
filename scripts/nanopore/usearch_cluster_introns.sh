#!/bin/bash

#SBATCH --partition sched_mit_cburge
#SBATCH --job-name=usearch_cluster_introns
#SBATCH --output=/home/ejkk/Cluster/IME/outfiles/usearch_cluster_filtered_introns_04.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=128000

echo "hola"

echo "clustering filtered non-concatemer introns"

usearch -cluster_fast /pool001/ejkk/RIL2/filtered_introns.fa -id 0.4 -consout /pool001/ejkk/RIL2/f_cons_04.fa -msaout /pool001/ejkk/RIL2/MSAs/f_msa_04.fa -uc /pool001/ejkk/RIL2/f_clusters_04.uc

echo "done!"

echo "merging MSAs into one file..."

tail -n +1 /pool001/ejkk/RIL2/MSAs/f_msa_04.fa* > /pool001/ejkk/RIL2/MSAs/f_msa_04_all.fa

echo "deleting individual msa files..."

rm /pool001/ejkk/RIL2/MSAs/f_msa_04.fa*

echo "all done :)"
