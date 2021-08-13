'''
Perform RNA-seq data QC using fastqc.
Quantify transcript expression in paired-end RNA-seq data with salmon.
Determine differential expression with voom.
Execute rna-seq_analysis.ipynb to summarize results.

Requirements: 
- update config.yaml with your paths and with conditions for each sample
- conda packages: fastqc, salmon, snakemake
- R bioconductor packages: tximport, ensembldb, EnsDb.Hsapiens.v79, edgeR, limma
- FASTA file containing your reference transcript
- FASTQ file(s) containing your reads

Usage: snakemake -j
Include -j to run with number of available CPU cores in the machine

dry run:
snakemake -np

work flow diagram:
snakemake --forceall --dag | dot -Tpng > dag.png
'''

import json
from os.path import join, basename, dirname
import pandas as pd
import re

# globals ----------------------------

configfile: 'config.yaml'

# Full path to an uncompressed FASTA file of reference transcriptome
CDNA = config['CDNA']

# Full path to data directory with fastq files
DATA_DIR = config['DATA_DIR']

# Full path to a folder where intermediate output files will be created.
OUT_DIR = config['OUT_DIR']

# conditions for the samples
CONDS = config['CONDS']

# path to trimmomatic installation
TRIM = config['TRIM']


# FUNCTION: sample_dictionary
# PARAMETER: fastq -> path to a folder with fastq files (FASTQ_DIR)
# RETURN: 
    # dictionary -> key: sample name ; value: dictionary where keys are 'R1' and/or 'R2' and values are the sample's full path
    # list_of_R1_and_R2 --> list of all samples without '_001.fastq.gz' ending

dictionary = {}
list_of_R1_and_R2 = []

def sample_dictionary(fastq):
	for root, dirs, files in os.walk(fastq):
		global dictionary
		global list_of_R1_and_R2
		for sample_files in files:
			if bool(re.search('_R._001.fastq.gz', sample_files)): # added because directory contains I1 and I1 files in addition to R1, R2
				sample_complete = sample_files.split('_001.fastq.gz')
				sample_info = sample_files.split('_')
				sample_name = sample_info[0] + '_' + sample_info[1]
				if sample_name in dictionary: 
					full_path = root + '/' + sample_files
					dictionary[sample_name][sample_info[4]] = [full_path]
					list_of_R1_and_R2.append(full_path)
				else: 
					dictionary[sample_name] = {}
					full_path = root + '/' + sample_files
					dictionary[sample_name][sample_info[4]] = [full_path]
					list_of_R1_and_R2.append(full_path) 
		return dictionary, list_of_R1_and_R2

FILES, SAMPLES_FULL_PATH = sample_dictionary(DATA_DIR)
SAMPLES = sorted(FILES.keys())

sample_names = []

def R1_R2_sample_names(list_r1_r2):
    global sample_names
    for i in range(0, len(list_r1_r2)):
        name_full = list_r1_r2[i].split('/')[-1]
        name = name_full.rstrip('.fastq.gz')
        sample_names.append(name)
    return sorted(sample_names)

SAMPLES_FULL = R1_R2_sample_names(SAMPLES_FULL_PATH)

# make sampleKey.csv needed for differential expression analysis
sampleKey = pd.DataFrame(
	{
		'sample': SAMPLES,
		'cond': CONDS
	}
)
sampleKey.to_csv('sampleKey.csv')

# Rules ------------------------------

rule all:
	input:
		join(dirname(CDNA), 'salmon', basename(CDNA).rstrip(".fa")),
		[OUT_DIR + "/" + x for x in expand('fastQC_output/{sample_full}_fastqc.html', sample_full = SAMPLES_FULL)],

		# trimmomatic output
		[OUT_DIR + "/" + x for x in expand('trimmedFQ/{sample}_filtered_1P.fastq.gz', sample = SAMPLES)],
		[OUT_DIR + "/" + x for x in expand('trimmedFQ/{sample}_filtered_1U.fastq.gz', sample = SAMPLES)],
		[OUT_DIR + "/" + x for x in expand('trimmedFQ/{sample}_filtered_2P.fastq.gz', sample = SAMPLES)],
		[OUT_DIR + "/" + x for x in expand('trimmedFQ/{sample}_filtered_2U.fastq.gz', sample = SAMPLES)],

		# 2nd round of fastqc
		[OUT_DIR + "/" + x for x in expand('fastQC_output_2/{sample}_filtered_1P_fastqc.html', sample = SAMPLES)],
		[OUT_DIR + "/" + x for x in expand('fastQC_output_2/{sample}_filtered_2P_fastqc.html', sample = SAMPLES)],

		[OUT_DIR + "/" + x for x in expand('{sample}/quant.sf', sample = SAMPLES)],
		[OUT_DIR + "/" + x for x in expand('{sample}/lib_format_counts.json', sample = SAMPLES)],
		OUT_DIR + '/quant.tsv.gz',

		OUT_DIR + '/counts_by_gene.csv',
		OUT_DIR + '/dispersion.png',
		OUT_DIR + '/de_summary.csv',

		# multiqc
		OUT_DIR + '/quality_control_metrics/multiqc/multiqc_report.html'

rule fast_qc_1:
	input:
		lambda wildcards: [s for s in SAMPLES_FULL_PATH if wildcards.sample_full in s]
	output:
		'{OUT_DIR}/fastQC_output/{sample_full}_fastqc.html'
	log:
		'{OUT_DIR}/logs/{sample_full}_fastqc.log',
	shell:
		"""
		fastqc -o {OUT_DIR}/fastQC_output {input} &> {log}
		"""

# trim fastq's to remove illumina adaptors
rule trim:
	input:
		r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
		r2 = lambda wildcards: FILES[wildcards.sample]['R2']
	output:
		'{OUT_DIR}/trimmedFQ/{sample}_filtered_1P.fastq.gz',
		'{OUT_DIR}/trimmedFQ/{sample}_filtered_1U.fastq.gz',
		'{OUT_DIR}/trimmedFQ/{sample}_filtered_2P.fastq.gz',
		'{OUT_DIR}/trimmedFQ/{sample}_filtered_2U.fastq.gz'
	log:
		summary = '{OUT_DIR}/logs/{sample}_summary_trimmomatic.log',
		detail = '{OUT_DIR}/logs/{sample}_detail_trimmomatic.log'
	shell:
		"""
		java -jar {TRIM}/trimmomatic-0.36.jar PE -trimlog {log.detail} {input.r1} {input.r2} -baseout {OUT_DIR}/trimmedFQ/{wildcards.sample}_filtered.fastq.gz ILLUMINACLIP:{TRIM}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &> {log.summary}
		"""

# run trimmed read 1's through fastQC
rule fast_qc_2_a:
	input:
		expand(join(OUT_DIR, 'trimmedFQ', '{sample}_filtered_1P.fastq.gz'), sample= SAMPLES)
	output:
		'{OUT_DIR}/fastQC_output_2/{sample}_filtered_1P_fastqc.html'
	log:
		'{OUT_DIR}/logs/{sample}_filtered_1P_fastqc_2.log'
	shell:
		"""
		fastqc -o {OUT_DIR}/fastQC_output_2 {OUT_DIR}/trimmedFQ/{wildcards.sample}_filtered_1P.fastq.gz &> {log}
		"""

# run trimmed read 2's through fastQC
rule fast_qc_2_b:
	input:
		expand(join(OUT_DIR, 'trimmedFQ', '{sample}_filtered_2P.fastq.gz'), sample= SAMPLES)
	output:
		'{OUT_DIR}/fastQC_output_2/{sample}_filtered_2P_fastqc.html'
	log:
		'{OUT_DIR}/logs/{sample}_filtered_2P_fastqc_2.log'
	shell:
		"""
		fastqc -o {OUT_DIR}/fastQC_output_2 {OUT_DIR}/trimmedFQ/{wildcards.sample}_filtered_2P.fastq.gz &> {log}
		"""

rule salmon_index:
	input:
		cdna = CDNA
	output:
		index = join(dirname(CDNA), 'salmon', basename(CDNA).rstrip(".fa"))
	log:
		 'logs/salmon_index.log'
	threads: 4
	shell:
		"""	
		salmon index -p {threads} -t {input} -i {output} &> {log}
		"""

rule salmon_quant:
	input:
		r1 = expand(join(OUT_DIR, 'trimmedFQ', '{sample}_filtered_1P.fastq.gz'), sample= SAMPLES),
		r2 = expand(join(OUT_DIR, 'trimmedFQ', '{sample}_filtered_2P.fastq.gz'), sample= SAMPLES),
		index = rules.salmon_index.output.index
	output:
		join(OUT_DIR, '{sample}', 'quant.sf'),
		join(OUT_DIR, '{sample}', 'lib_format_counts.json')#,
		#join(OUT_DIR, '{sample}', 'map.sam')
	log:
		'logs/{sample}_salmons_quant.log'
	threads: 4
	shell:
		"""
		salmon quant -p {threads} --validateMappings -i {input.index} -l A -1 <(gunzip -c {OUT_DIR}/trimmedFQ/{wildcards.sample}_filtered_1P.fastq.gz) -2 <(gunzip -c {OUT_DIR}/trimmedFQ/{wildcards.sample}_filtered_2P.fastq.gz) -o {OUT_DIR}/{wildcards.sample} &> {log}
		"""

rule collate_salmon:
    input:
        expand(join(OUT_DIR, '{sample}', 'quant.sf'), sample= SAMPLES)
    output:
        '{OUT_DIR}/quant.tsv.gz'
    run:
        import gzip

        b = lambda x: bytes(x, 'UTF8')

        # Create the output file.
        with gzip.open(output[0], 'wb') as out:

            # Print the header.
            header = open(input[0]).readline()
            out.write(b('sample\t' + header))

            for i in input:
                sample = basename(dirname(i))
                lines = open(i)
                # Skip the header in each file.
                lines.readline()
                for line in lines:
                    out.write(b(sample + '\t' + line))


# run MultiQC
rule multi_qc_metrics:
    input:
        fastqc_files = expand(join(OUT_DIR, 'fastQC_output', '{sample_full}_fastqc.html'), sample_full = SAMPLES_FULL),
        salmon_files = expand(join(OUT_DIR, '{sample}', 'quant.sf'), sample = SAMPLES) 
    output:
        '{OUT_DIR}/quality_control_metrics/multiqc/multiqc_report.html'
    shell:
    	"""
		multiqc {OUT_DIR} -o '{OUT_DIR}/quality_control_metrics/multiqc'
		"""

#differential expression analysis
rule diff_expression:
	input:
		expand(join(OUT_DIR, '{sample}', 'quant.sf'), sample= SAMPLES),
		'config.yaml'
	output:
		join(OUT_DIR, 'counts_by_gene.csv'),
		join(OUT_DIR, 'dispersion.png'),
		join(OUT_DIR, 'de_summary.csv'),
		join(OUT_DIR, 'diff_all.csv')
	shell:
		"""
		Rscript --vanilla scripts/de.R {OUT_DIR} sampleKey.csv
		"""