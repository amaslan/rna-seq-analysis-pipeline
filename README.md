# rna-seq-pipeline
## Website: https://amaslan.github.io/rna-seq-pipeline-ivana/

## Pipeline summary
1. FASTQC - high-level sequencing quality control
2. TRIMMOMATIC - remove adapters
3. FASTQC - round 2
4. SALMON - pseudo-mapping & transcript quantification
5. PICARD - insert size distribution, GC bias, alignment summary
6. MULTIQC
7. VOOM/LIMMA - differential expression analysis

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
    $ conda install nbconvert
    $ conda install ipykernel
    $ conda install scikit-learn
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
module load ...
screen -S <job_name> snakemake -j
```
## Key output

### Sequencing quality control metics & plots
- fastqc.html files: These are important for determining whether pre-processing of the reads is needed before differential expression analysis. Key statistics analyzed are per base sequence quality, per tile sequence quality, per sequence quality scores, per base sequence content, per sequence GC content, per base N content, sequence length distribution, sequence duplication levels, overrepresented sequences, and adapter content.

### Picard metrics & plots
- multiple_metrics.alignment_summary_metrics: high-level summary of alignment statistics like total PF reads and percent adapter.
- multiple_metrics.gc_bias.detail_metrics: data used to generate coverage vs. GC content curves in rna-seq_analysis.ipynb.
- multiple_metrics.gc_bias.pdf: normalized coverage vs. GC content with distribution for samples and reference transcriptome indicated.
- multiple_metrics.gc_bias.summary_metrics: statistics describing GC bias such as GC dropout.
- multiple_metrics.insert_size_histogram.pdf: histogram of insert sizes
- multiple_metrics.insert_size_metrics: data used to generate insert size histograms in rna-seq_analysis.ipynb.

### BAM file and transcript counts
- mapped.sorted.cleaned.bam: final, cleaned, aligned, sorted sequence file
- quant.tsv.gz: transcript quantification from Salmon structured as sample, transcript name, length, effective length, transcripts per million (TPM), and number of reads.

### Differential expression
- counts_by_gene.csv - gene-level counts from aggregating Salmon transcript counts.
- dispersion.png - output from limma's voom function showing noise model.
- de_summary.csv - summary of number genes that are more lowly expressed, similarly expressed, or more highly expressed for all pairwise comparisons.
- venn_de_genes.png - venn diagram showing the overlapping number of differentially expressed genes for all pairwise comparisons.
- diff_all.csv - gene id's for differentially expressed genes at the intersection of all pairwise comparisons. Use this gene list to perform GO gene set enrichment analysis for example.
- volcano plots - probability that a gene is differentially expressed vs. the log fold change in expression.
- histograms of p-values - distribution of p-values for differential expression analysis (ensure not uniform at peak for p < 0.05).

### Jupyter notebook plots
- isize_normalized.png - plot of each sample's insert size distribution
- gc_bias.png - plot of each sample's coverage as a funciton of GC content
- uniquetx_vs_reads.png - count of unique transcripts vs. total reads for each sample  
- de_heatmap.png - heatmap of counts by gene for genes that are differentially expressed for all pairwise comparisons
- pca.png - principal component analysis using counts by gene for genes that are differentially expressed for all pairwise comparisons

## Known issues
- rna-seq_metrics.ipynb is intended as a model for looking at the sequencing metrics output from the sequencer. Because the output format varies by sequencer, this notebook must be adjusted. However, the metrics examined are important to consider for all sequencing runs; the initial data formatting and input needs to be updated by each user.
- rna-seq_analysis.ipynb aggregates all of the pipeline results. Embedding images yields issues when converting to HTML. All embedded images can be seen and everything works well when running the notebook. However, to save as HTML with the embedded images, the current workaround is to use 'File --> Print Preview.' This is a known bug with image embedding and I'm working to find a better fix.

## References
In creating the Snakefile for this pipeline, the below two repositories were very helpful for modeling the pipeline framework.
- https://github.com/crazyhottommy/RNA-seq-analysis/tree/master/RNA-seq-snakemake-pipeline
- https://github.com/slowkow/snakefiles