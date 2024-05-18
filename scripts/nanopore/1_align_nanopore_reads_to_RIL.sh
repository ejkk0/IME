#!/bin/bash

#SBATCH --partition sched_mit_cburge
#SBATCH --job-name=align_nanopore_reads_to_RIL_ORI_sam
#SBATCH --output=/nobackup1/ejkk/Cluster/IME/RIL2/outfiles/minimap_sam.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=128000

./minimap2/minimap2 -ax map-ont /nobackup1/ejkk/Cluster/IME/RIL2/RIL_genome_ORI_start.fa /nobackup1/ejkk/Cluster/IME/RIL2/all_nanopore_reads.fastq.gz > /nobackup1/ejkk/Cluster/IME/RIL2/RIL_ORI_alignment.sam
