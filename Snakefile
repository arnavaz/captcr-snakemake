#### This Snakemake file is to generare mixcr clones files from the fastqs###
snakefile: "/path/to/extract.smk"
configfile: "config.yaml"
###values for running mixcr### these should be set in the config file###

java_align = "java -jar  mixcr.jar  align -p rna-seq -s hsa -OallowPartialAlignments=true -OvParameters.geneFeatureToAlign=VGeneWithP -r "
assemble1 = "java -jar   mixcr.jar assemblePartial " 
assembleEx="java -jar   mixcr.jar extendAlignments "
export="java -jar  mixcr.jar exportClones"

aligndir=config['aligndir']
basedir= config['basedir']
clonetrack=config['clonetrack']
logdir=config['logdir']
###  Rules #####

rule all:
         input= expand("clone_track" + "{sample}.{param}.output.pdf", param=config["patterns"])



rule align:

	input:
	      R1 = expand("basedir/{id}_barcode_R1.fastq", id=IDS1),
	      R2 = expand("basedir/{id}_barcode_R1.fastq", id=IDS2)


	
	output: 
		vdjca = "aligndir/alignments_{id}.vdjca"
	shell:
		'''java_align +  logdir + "/{id}.txt " + basedir + "/{input.R1}" + basedir "/{input.R2}" + output.vdjca'''


rule assemble1:

	input:
	       vdjca = aligndir + "alignments_" + "{id}" + ".vdjca"

	output:
	       vdjca_res1= assembledir + "rescue1_" +"{id}"+ ".vdjca"
		
	shell:
               '''assemble1 + aligndir + "alignments_" + "{id}" + ".vdjca " + {output.vdjca_res1}'''



rule assemble2:

	input:
	       vdjca_res1 = "rescue1_" +"{id}"+ ".vdjca"

	output:
	        vdjca_res2=assembledir + "rescue2_" + "{id}" + ".vdjca"

	shell:

	       '''assemble1 + assembledir +  "rescue1_" + "{id}" + ".vdjca " + {output.vdjca_res2}'''

rule assembleEx:

	input:
		vdjca_res2=assembledir + "rescue2_" + "{id}" + ".vdjca"

	output:
               vdjca_ext=assembledir + "extended_" + "{id}" + ".vdjca"
	shell:
	      '''assembleEx + assembledir + "rescue2_" + "{id}" + ".vdjca " + {output.vdjca_ext}'''


rule assemble:

		input:
		      vdjca_ext=assembledir + "extended_" + "{id}" + ".vdjca"
	        output:
	              clns=assembledir  + "{id}" + ".clns",
			log=logdir + "log_assemble_" + "{id}" +".txt "		
		
		    

		shell:
		      '''assemble + logdir + "log_assemble_" + "{id}" +".txt " +  assembledir + "extended_" + "{id}" +
		      ".vdjca " + {output.clns} '''


rule export_A:
	    input:
		   clns=assembledir  + "{id}" + ".clns"
		
	    output:
          	   clones_A=  clonesdir + "CLONES_TRA" + "{id}"	
		
	
	shell:
		
		''' export "exportClones --chains TRA " + assembledir + "{id}" + ".clns " + {output.clones_A} '''
rule export_B:
	
	    input:
		   clns=assembledir  + "{id}" + ".clns"
	    output:
          	   clones_B=  clonesdir + "CLONES_TRB" + "{id}"		
	
	    shell:
		
		''' export "exportClones --chains TRB " + assembledir + "{id}" + ".clns " + {output.clones_B} ''' 		

rule log_pars:
	
	input:
		log=logdir + "log_assemble_" + "{id}" +".txt "
	output: 
		a1=logdir + "/assemble_stats.csv",
		a2=logdir + "/align_stats.csv"
	shell:
		
		''' python parse_log.py + "-i" +  logdir + "-o" + logdir'''
		
	
rule run_QC_plotter:
	
	input:
		a1=logdir + "/assemble_stats.csv",
		a2=logdir + "/align_stats.csv"
		
	output:
		qc="project_QC.pdf"
		
	shell:
		
		''' Rscript mixcr_QCplot.R {input.a1} {input.a2} projectname '''
rule run_clone_tracker_TRA:
	
	input:
		clones_A=clonesdir + "CLONES_TRA" + "{id}",
		clones_B=clonesdir + "CLONES_TRB" + "{id}",
	output:
		clone_tra = "clonetrack_" + pattern + {input.clones_A} + ".pdf",
		clone_trb = "clonetrack_" + pattern + {input.clones_B} + ".pdf",
		
		
		
	
	
	
