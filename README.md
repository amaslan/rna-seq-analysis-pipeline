# rna-seq-pipeline
## Website: https://amaslan.github.io/rna-seq-pipeline-ivana/

## Pipeline summary
1. FASTQC - high-level sequencing quality control
2. TRIMMOMATIC - remove adapters
3. FASTQC - round 2 after trimming
4. SALMON - pseudo-mapping & transcript quantification
5. MULTIQC
6. VOOM/LIMMA - differential expression analysis

## Steps
### 1. Input files required:
- FASTA file containing your reference transcripts - wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz --> from: http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/
- FASTA/FASTQ file(s) containing your reads

### 2. Installations:
- miniconda and bioconda
```
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
$ conda config --add channels bioconda
```

- snakemake, salmon, nbconvert, ipykernel, scikit-learn, etc.
    - Set up conda environment
    ```
    conda env create rna-seq-env
    source activate rna-seq-env
    ```
    - Key ones to install in your own environment:
    ```
    $ conda install -c bioconda snakemake
    $ conda install salmon
    ```

- R bioconductor packages for differential expression script. Installation via bioconductor depends on R version
```
# older
> source("https://bioconductor.org/biocLite.R")
> biocLite("tximport")
> biocLite("ensembldb")
> biocLite("edgeR")
> biocLite("limma")

# newer R 3.5 and up
> if (!requireNamespace("BiocManager", quietly = TRUE))
>     install.packages("BiocManager")

> BiocManager::install("tximport")
> BiocManager::install("ensembldb") #had to do outside conda env
> BiocManager::install("edgeR")
```

### 3. Update config.yaml with file paths and conditions
Also get and update gtf file for DE analysis! I stored here: (see de.R script)
"/clusterfs/rosalind/groups/streetslab/amaslan/ivana/annotations/Homo_sapiens.GRCh38.104.gtf.gz"

### 4. Run with the number of available CPU cores in the machine:
```
module load ... # whatever modules you need to load to run the pipeline
screen -S <job_name> snakemake -j # or submit a batch job with shell script like .sh file examples in repo
```
## Key output

### Sequencing quality control metics & plots
- multiqc.html file: Lets you see all the fastqc and salmon results for all samples in one place so you don't have to go through files individually.
- fastqc.html files: These are important for determining whether pre-processing of the reads is needed before differential expression analysis. Key statistics analyzed are per base sequence quality, per tile sequence quality, per sequence quality scores, per base sequence content, per sequence GC content, per base N content, sequence length distribution, sequence duplication levels, overrepresented sequences, and adapter content.

### transcript counts
- quant.tsv.gz: transcript quantification from Salmon structured as sample, transcript name, length, effective length, transcripts per million (TPM), and number of reads.

### Differential expression
- counts_by_gene.csv - gene-level counts from aggregating Salmon transcript counts.
- dispersion.png - output from limma's voom function showing noise model.
- de_summary.csv - summary of number genes that are more lowly expressed, similarly expressed, or more highly expressed for all pairwise comparisons.
- volcano plots - probability that a gene is differentially expressed vs. the log fold change in expression.
- histograms of p-values - distribution of p-values for differential expression analysis (ensure not uniform at peak for p < 0.05).

## References
In creating the Snakefile for this pipeline, the below two repositories were very helpful for modeling the pipeline framework.
- https://github.com/crazyhottommy/RNA-seq-analysis/tree/master/RNA-seq-snakemake-pipeline
- https://github.com/slowkow/snakefiles