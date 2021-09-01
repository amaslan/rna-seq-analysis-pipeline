#!/bin/bash
# Job name:
#SBATCH --job-name=amm_rna
#
# Account:
#SBATCH --account=co_rosalind
#
# Partition:
#SBATCH --partition=savio2_htc
#
# Quality of Service:
#SBATCH --qos=rosalind_htc2_normal
#
# Wall clock limit:
#SBATCH --time=1-00:00:00
#
# Nodes:
#SBATCH --nodes=1
# Specify number of tasks for use case (example):
#SBATCH --ntasks-per-node=10
## Command(s) to run:

cd /clusterfs/rosalind/groups/streetslab/amaslan/ivana/data_ref


wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR590/ERR590408/ERR590408_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR590/ERR590408/ERR590408_2.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR590/ERR590410/ERR590410_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR590/ERR590410/ERR590410_2.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR590/ERR590401/ERR590401_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR590/ERR590401/ERR590401_2.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR590/ERR590400/ERR590400_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR590/ERR590400/ERR590400_2.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR590/ERR590398/ERR590398_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR590/ERR590398/ERR590398_2.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR590/ERR590399/ERR590399_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR590/ERR590399/ERR590399_2.fastq.gz