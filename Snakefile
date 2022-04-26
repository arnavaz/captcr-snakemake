#### This Snakemake file is to generare mixcr clones files from the fastqs###

configfile: "path/to/config.yaml"
###values for running mixcr### these should be set in the config file###

java_align = "java -jar $mixcr_dir/mixcr.jar  align -p rna-seq -s hsa -OallowPartialAlignments=true -OvParameters.geneFeatureToAlign=VGeneWithP -r "
assemble1 = "java -jar $mixcr_dir/mixcr.jar assemblePartial " 
assembleEx="java -jar $mixcr_dir/mixcr.jar extendAlignments "
export="java -jar $mixcr_dir/mixcr.jar exportClones"
aligndir=

###  Rules #####

rule all:
         input= expand("clone_track" + "{sample}.{param}.output.pdf", sample=config["samples"], param=config["patterns"])



rule align:

	input:
		fastq1=cwd + "/{id}_barcode_R1.fastq",
		fastq2=cwd + "/{id}_barcode_R2.fastq",


	
	output: 
		vdjca = aligndir + "alignments_" + "{id}" + ".vdjca"
	shell:
		'''java_align +  logdir + "/{id}.txt " + cwd  + "/{input.fastq1}" + cwd "/{input.fastq2}" + output.vdjca	'''


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
	              clns=assembledir  + "{id}" + ".clns"

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
	
