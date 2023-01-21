rule run_ASCAT:
	params:
		threads="8",
		runtime="4:00",
		memory="16000",
		sex = lambda wcs: df_samplesheet.loc[wcs.sample,"sex"]
	input:
		BAM= config["working_dir"]+"/BAM/{sample}.bam"
	output:
		ratios=config["working_dir"]+"/out_nocontrol/CNA/ASCAT/{sample}/{sample}_segments.tsv",
		CNAs=config["working_dir"]+"/out_nocontrl/CNA/ASCAT/{sample}/{sample}.ASCATprofile.png"
	conda:
		"../envs/WGS.yaml"
	shell:
		"Rscript scripts/ezASCAT.R {input.BAM} {wildcards.sample} {config[working_dir]}/out_nocontrol/CNA/ASCAT/{wildcards.sample};"\
		"mv {config[working_dir]}/out_nocontrol/CNA/ASCAT/{wildcards.sample}/logR.ASCATprofile.png {config[working_dir]}/out_nocontrol/CNA/ASCAT/{wildcards.sample}/{wildcards.sample}.ASCATprofile.png;"
		