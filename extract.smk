
import os
import sys
configfile: "config.yaml"
basedir =
IDS1, = glob_wildcards(input{id}.fastq.gz")
IDS2, = glob_wildcards("thedir/{id}.fastq.gz")


rule all:
    input:  expand(basedir + "/" +{id}_barcode_R1.fastq", id=IDS1),
	    expand(basedir + "/"+{id}_barcode_R2.fastq", id=IDS2)


rule extract:

        conda: 'conda_env/biopython.yml'
	input: 
              R1=otherdir/{id}_R1.fastq, 
	      R2=otherdir/{id}_R2.fastq
        output:
	     R1 = cwd + "{id}_barcode_R1.fastq",
	     R2 = cwd + "{id}_barcode_R2.fastq"


       outprefix = lambda wildcards, output: output.R1.split('_barcode_')[0]
       shell:
        '''
        python {extract_barcodes}
            --read1 {{input.R1}}
            --read2 {{input.R2}}
            --outfile {{params.outprefix}}
            {{params.barcodes}}
        '''.format(extract_barcodes = config['paths']['dependencies']['extract_barcodes_path']

