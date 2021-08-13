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

- R bioconductor packages for differential expression script
```
> source("https://bioconductor.org/biocLite.R")
> biocLite("tximport")
> biocLite("ensembldb")
> biocLite("edgeR")
> biocLite("limma")
```

### 3. Update config.yaml with file paths and conditions

### 4. Run with the number of available CPU cores in the machine:
```
module load ... # whatever modules you need to load to run the pipeline
screen -S <job_name> snakemake -j
```
## Key output

### Sequencing quality control metics & plots
- fastqc.html files: These are important for determining whether pre-processing of the reads is needed before differential expression analysis. Key statistics analyzed are per base sequence quality, per tile sequence quality, per sequence quality scores, per base sequence content, per sequence GC content, per base N content, sequence length distribution, sequence duplication levels, overrepresented sequences, and adapter content.

### BAM file and transcript counts
- mapped.sorted.cleaned.bam: final, cleaned, aligned, sorted sequence file
- quant.tsv.gz: transcript quantification from Salmon structured as sample, transcript name, length, effective length, transcripts per million (TPM), and number of reads.

### Differential expression
- counts_by_gene.csv - gene-level counts from aggregating Salmon transcript counts.
- dispersion.png - output from limma's voom function showing noise model.
- de_summary.csv - summary of number genes that are more lowly expressed, similarly expressed, or more highly expressed for all pairwise comparisons.
- diff_all.csv - gene id's for differentially expressed genes at the intersection of all pairwise comparisons. Use this gene list to perform GO gene set enrichment analysis for example.
- volcano plots - probability that a gene is differentially expressed vs. the log fold change in expression.
- histograms of p-values - distribution of p-values for differential expression analysis (ensure not uniform at peak for p < 0.05).

## References
In creating the Snakefile for this pipeline, the below two repositories were very helpful for modeling the pipeline framework.
- https://github.com/crazyhottommy/RNA-seq-analysis/tree/master/RNA-seq-snakemake-pipeline
- https://github.com/slowkow/snakefiles