#!/bin/bash

#SBATCH --partition sched_mit_cburge
#SBATCH --job-name=usearch_cluster_barcodes
#SBATCH --output=/nobackup1/ejkk/Cluster/IME/RIL2/usearch_cluster_barcodes.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=128000

usearch -cluster_smallmem /nobackup1/ejkk/Cluster/IME/RIL2/nanopore_full_barcode_list.fa -id 0.8 -sortedby other -clusters /nobackup1/ejkk/Cluster/IME/RIL2/barcode_clusters/c_