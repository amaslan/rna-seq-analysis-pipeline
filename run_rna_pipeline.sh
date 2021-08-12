#!/bin/bash
# Job name:
#SBATCH --job-name=amm_rna_dam
#
# Account:
#SBATCH --account=co_rosalind
#
# Partition:
#SBATCH --partition=savio
#
# Quality of Service:
#SBATCH --qos=rosalind_savio_normal
#
# Wall clock limit:
#SBATCH --time=3-00:00:00
#
# Nodes:
#SBATCH --nodes=1
## Command(s) to run:
module load fastqc
module load samtools/1.8
module load multiqc
snakemake -j --rerun-incomplete