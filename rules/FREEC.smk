rule run_FREEC:
	params:
		threads="8",
		runtime="4:00",
		memory="16000",
		sex = lambda wcs: df_samplesheet.loc[wcs.sample,"sex"]
	input:
		BAM= config["working_dir"]+"/BAM/{sample}.bam"
	output:
		ratios=config["working_dir"]+"/out_nocontrol/CNA/FREEC/{sample}/{sample}.bam_ratio.txt",
		CNAs=config["working_dir"]+"/out_nocontrol/CNA/FREEC/{sample}/{sample}.bam_CNVs"
	conda:
		"../envs/WGS.yaml"
	shell:
		"cp {config[template_FREEC_nocontrol]} {config[working_dir]}/out_nocontrol/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{chrLenFile}}#{config[chrLenFile]}#g' {config[working_dir]}/out_nocontrol/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{outdir}}#{config[working_dir]}/out_nocontrol/CNA/FREEC/{wildcards.sample}#g' {config[working_dir]}/out_nocontrol/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{BAMFILE}}#{input.BAM}#g' {config[working_dir]}/out_nocontrol/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{SEX}}#{params.sex}#g' {config[working_dir]}/out_nocontrol/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"{config[freec_bin]} -config {config[working_dir]}/out_nocontrol/CNA/FREEC/{wildcards.sample}/config.txt"


rule run_FREEC_BAF:
	params:
		threads="8",
		runtime="4:00",
		memory="16000",
		sex = lambda wcs: df_samplesheet.loc[wcs.sample,"sex"]
	input:
		BAM= config["working_dir"]+"/BAM/{sample}.bam"
	output:
		ratios=config["working_dir"]+"/out_nocontrol/CNA/FREEC_BAF/{sample}/{sample}.bam_ratio.txt",
		CNAs=config["working_dir"]+"/out_ncontrol/CNA/FREEC_BAF/{sample}/{sample}.bam_CNVs"
	conda:
		"../envs/WGS.yaml"
	shell:
		"cp {config[template_FREEC_BAF_nocontrol]} {config[working_dir]}/out_nocontrol/CNA/FREEC_BAF/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{chrLenFile}}#{config[chrLenFile]}#g' {config[working_dir]}/out_nocontrol/CNA/FREEC_BAF/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{outdir}}#{config[working_dir]}/out_nocontrol/CNA/FREEC_BAF/{wildcards.sample}#g' {config[working_dir]}/out_nocontrol/CNA/FREEC_BAF/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{BAMFILE}}#{input.BAM}#g' {config[working_dir]}/out_nocontrol/CNA/FREEC_BAF/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{SEX}}#{params.sex}#g' {config[working_dir]}/out_nocontrol/CNA/FREEC_BAF/{wildcards.sample}/config.txt;"\
		"{config[freec_bin]} -config {config[working_dir]}/out_nocontrol/CNA/FREEC_BAF/{wildcards.sample}/config.txt"

rule run_FREEC_control:
	params:
		threads="8",
		runtime="4:00",
		memory="16000",
		sex = lambda wcs: df_samplesheet.loc[wcs.sample,"sex"]
	input:
		BAM_normal= config["working_dir"]+"/BAM/{sample}_control.bam",
		BAM_tumor= config["working_dir"]+"/BAM/{sample}.bam"
	output:
		ratios=config["working_dir"]+"/out_control/CNA/FREEC/{sample}/{sample}.bam_ratio.txt",
		CNAs=config["working_dir"]+"/out_control/CNA/FREEC/{sample}/{sample}.bam_CNVs"
	conda:
		"../envs/WGS.yaml"
	shell:
		"cp {config[template_FREEC_control]} {config[working_dir]}/out_control/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{chrLenFile}}#{config[chrLenFile]}#g' {config[working_dir]}/out_control/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{outdir}}#{config[working_dir]}/out_control/CNA/FREEC/{wildcards.sample}#g' {config[working_dir]}/out_control/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{BAMFILE_NORMAL}}#{input.BAM_normal}#g' {config[working_dir]}/out_control/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{BAMFILE_TUMOR}}#{input.BAM_tumor}#g' {config[working_dir]}/out_control/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"sed -i 's#{{SEX}}#{params.sex}#g' {config[working_dir]}/out_control/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"{config[freec_bin]} -config {config[working_dir]}/out_control/CNA/FREEC/{wildcards.sample}/config.txt;"\
		"cp {config[working_dir]}/out_control/CNA/FREEC/{wildcards.sample}/{wildcards.sample}.bam_ratio.txt {output.ratios}; "\
		"cp {config[working_dir]}/out_control/CNA/FREEC/{wildcards.sample}/{wildcards.sample}.bam_CNVs {output.CNAs}; "\

rule create_vizualization_datafiles_FREEC:
	params:
		threads="1",
		runtime="40",
		memory="8000"
	input:
		ratios = config["working_dir"]+"/{out}/CNA/FREEC/{sample}/{sample}.bam_ratio.txt",
		CNAs = config["working_dir"]+"/{out}/CNA/FREEC/{sample}/{sample}.bam_CNVs",
	output:
		circos = config["working_dir"]+"/{out}/CNA/FREEC/{sample}/{sample}_CNA_circos.txt",
		chrplot = config["working_dir"]+"/{out}/CNA/FREEC/{sample}/{sample}_CNA_chrplots.tsv"
	conda:
		"../envs/WGS.yaml"
	shell:
		"python scripts/freec2circos.py -i {input.ratios} -o {output.circos} --cnv {input.CNAs};"\
		"python scripts/freec2chrplots.py -i {input.ratios} -o {output.chrplot} --cnv {input.CNAs};"

rule FREEC_to_shatterseek:
	params:
		threads="1",
		runtime="20",
		memory="8000"
	input:
		config["working_dir"]+"/{out}/CNA/FREEC/{sample}/{sample}.bam_ratio.txt"
	output:
		config["working_dir"]+"/{out}/CNA/FREEC/{sample}/{sample}_shatterseekCNA.tsv"

	conda:
		"../envs/WGS.yaml"
	shell:
		"python scripts/freec2shatterseek.py -i {input} -o {output}"