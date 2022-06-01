#### This Snakemake file is to generate mixcr clones tracking figures from the fastqs###
configfile: "config.yaml"
basedir = config['basedir']
inputdir = config['inputdir']
aligndir=config['aligndir']
clonetrack=config['clonetrack']
logdir=config['logdir']
pars_log=config['pars_log_path']
mixcr_path=config['mixcr_path']
###values for running mixcr### these should be set in the config file###

java_align = "java -jar  mixcr.jar  align -p rna-seq -s hsa -OallowPartialAlignments=true -OvParameters.geneFeatureToAlign=VGeneWithP -r "
assemble1 = "java -jar   mixcr.jar assemblePartial " 
assembleEx="java -jar   mixcr.jar extendAlignments "
assemble_f="java -jar /cluster/tools/software/centos7/mixcr/3.0.12/mixcr.jar assemble -r"
export="java -jar  mixcr.jar exportClones"

IDS1, = glob_wildcards("inputdir/{id}_R1.fastq.gz")
IDS2, = glob_wildcards("inputdir/{id}_R2.fastq.gz")

###  Rules ####

rule all:
    input: 
	  expand("outputdir/clone_track_{id}.{param}.pdf", param=config["patterns"])

		

rule extract:
	
        conda: 'conda_env/biopython.yml'
	input: 
              R1=expand("inputdir/{id}_R1.fastq.gz", id=IDS1), 
	      R2=expand("inputdir/{id}_R2.fastq.gz",id=IDS2)
        output:
	     R1 = expand("basedir/{id}_barcode_R1.fastq", id=IDS1),
	     R2 = expand("basedir/{id}_barcode_R1.fastq", id=IDS2)

       params:
            outprefix = lambda wildcards, output: output.R1.split('_barcode_')[0]
            barcodes = config['barcodes']
       shell:
        '''
        python {extract_barcodes}
            --read1 {{input.R1}}
            --read2 {{input.R2}}
            --outfile {{params.outprefix}}
            {{params.barcodes}}
        '''.format(extract_barcodes = config['extract_barcodes_path'])


rule align:

	input:
	      R1 = expand("basedir/{id}_barcode_R1.fastq", id=IDS1),
	      R2 = expand("basedir/{id}_barcode_R1.fastq", id=IDS2)


	
	output: 
		vdjca = "aligndir/alignments_{id}.vdjca"
	shell:
		'''{java}  
		   logdir/{id}.txt  
		   basedir/{input.R1}
		   basedir/{input.R2} 
		   output.vdjca
		   '''.format(java=java_align)


rule assemble1:

	input:
	       vdjca = "aligndir/alignments_{id}.vdjca"

	output:
	       vdjca_res1= "assembledir/rescue1_{id}.vdjca"
		
	shell:
               '''{assemble}
	           aligndir/alignments_{id}.vdjca 
		  {output.vdjca_res1}
	       '''.format(assemble=assemble1)



rule assemble2:

	input:
	       vdjca_res1 = "rescue1_{id}.vdjca"

	output:
	        vdjca_res2=assembledir/rescue2_{id}.vdjca"

	shell:

	       '''{assemble}
	          assembledir/rescue1_{id}.vdjca 
		  {output.vdjca_res2}
	       '''.format(assemble=assemble1)

rule assembleEx:

	input:
		vdjca_res2=assembledir/rescue2_{id}.vdjca"

	output:
               vdjca_ext=assembledir/extended_{id}.vdjca"
	shell:
	      '''{assemble}
	         assembledir/rescue2_{id}.vdjca 
		 {output.vdjca_ext}
	      '''.format(assemble=assembleEx)


rule assemble:

	input:
		vdjca_ext=assembledir/extended_{id}.vdjca"
	output:
	        clns=assembledir/{id}.clns",
		log=logdir/log_assemble_id}.txt "			    

	shell:
	      '''{assemble}
		 logdir/log_assemble_{id}.txt 
	         assembledir/extended_{id}.vdjca
	         {output.clns} 
	     '''.format(assemble=assemble_f)


rule export_A:
	    input:
		   clns="assembledir/{id}.clns"
		
	    output:
          	   clones_A= "clonesdir/CLONES_TRA_{id}.txt"	
		
	
	shell:
		
		'''export 
		  exportClones --chains TRA
		  assembledir/{id}.clns  {output.clones_A}
		'''
rule export_B:
	
	    input:
		   clns="assembledir/{id}.clns"
	    output:
          	   clones_B="clonesdir/CLONES_TRB_{id}.txt"		
	
       shell:
		
		'''export 
		exportClones --chains TRB
		assembledir  {id}.clns
	        {output.clones_B}
	        ''' 		

rule log_pars:
	
	output: 
		a1="logdir/assemble_stats.csv",
		a2="logdir/align_stats.csv"
	shell:
		
		''' python {parse_log} 
		           -i  logdir  -o logdir
		'''.format(pars_log=confgi['pars_log_path'])
		
rule run_QC_plotter:
	
	input:
		a1="logdir/assemble_stats.csv",
		a2="logdir/align_stats.csv"
		
	output:
		qc="project_QC.pdf"
		
	shell:
		
		''' Rscript {mixcr_qc}
		    {input.a1} 
		    {input.a2} 
	        '''
rule  clone_tracker:
	
	input:
		clones='''clonesdir/CLONES_{chain}_{id}.txt'''.format(chain=config['chain'])
		
	output:
		clones_pdf="clonesdir/clone_track_{id}.{param}.pdf".format(param=config['patterns'])
		
	shell:
		'''Rscript {clonetracker} 
		   clonesdir
		   {param}
		   {type}
		'''.format(param=config['pattern'], type= config['chain])  
		   
		
		
	
	
	
