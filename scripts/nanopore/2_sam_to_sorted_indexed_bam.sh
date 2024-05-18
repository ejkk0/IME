#!/bin/bash

#SBATCH --partition sched_mit_cburge
#SBATCH --job-name=sam2bam
#SBATCH --output=/home/ejkk/Cluster/IME/outfiles/sort_bam.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=64000


echo "hihihi"
module load engaging/samtools/1.3.1
echo "converting sam -> bam..."
#samtools view -S -b $1.sam > $1.bam
samtools calmd -u /nobackup1/ejkk/Cluster/IME/RIL2/RIL_ORI_alignment.sam /nobackup1/ejkk/Cluster/IME/RIL2/RIL_genome_ORI_start.fa > /pool001/ejkk/RIL2/RIL_ORI_alignment_2.bam
echo "sorting bam..."
samtools sort /pool001/ejkk/RIL2/RIL_ORI_alignment_2.bam -o /pool001/ejkk/RIL2/RIL_ORI_alignment_2.sorted.bam
echo "indexing bam..."
samtools index /pool001/ejkk/RIL2/RIL_ORI_alignment_2.sorted.bam
echo "yay done"