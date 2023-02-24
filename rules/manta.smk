rule run_manta_nocontrol:
	params:
		threads="8",
		runtime="20:00",
		memory="8000"
	input:
		BAM= config["working_dir"]+"/BAM/{sample}.bam"
	output:
		config["working_dir"]+"/out_nocontrol/SV/manta/{sample}/{sample}_SV_unfiltered.vcf"
	conda:
		"../envs/manta.yaml"
	shell:
		"python {config[manta_bin]}/configManta.py --tumorBam {input.BAM} "\
		"--referenceFasta {config[reference_fasta]} "\
		"--runDir ${{TMPDIR}}/manta_nocontrol/{wildcards.sample} "\
		"--callRegions data/hs37d5_chr.bed.gz --exome;"\
		"python ${{TMPDIR}}/manta_nocontrol/{wildcards.sample}/runWorkflow.py -j {params.threads} -g 8; "\
		"cp ${{TMPDIR}}/manta_nocontrol/{wildcards.sample}/results/variants/tumorSV.vcf.gz {config[working_dir]}/out_nocontrol/SV/manta/{wildcards.sample}/{wildcards.sample}_SV_unfiltered.vcf.gz; "\
		"gunzip -d {config[working_dir]}/out_nocontrol/SV/manta/{wildcards.sample}/{wildcards.sample}_SV_unfiltered.vcf.gz; "\
		"rm -r ${{TMPDIR}}/manta_nocontrol/{wildcards.sample}"


#"--runDir {config[working_dir]}/tmp/SV/manta_nocontrol/{wildcards.sample} "\

rule run_manta_control:
	params:
		threads="8",
		runtime="20:00",
		memory="16000"
	input:
		BAM_normal= config["working_dir"]+"/BAM/{sample}_control.bam",
		BAM_tumor= config["working_dir"]+"/BAM/{sample}.bam",
	output:
		config["working_dir"]+"/out_control/SV/manta/{sample}/{sample}_SV_somatic.vcf"
	conda:
		"../envs/manta.yaml"
	shell:
		"python {config[manta_bin]}/configManta.py --normalBam {input.BAM_normal} --tumorBam {input.BAM_tumor} "\
		"--referenceFasta {config[reference_fasta]} "\
		"--runDir ${{TMPDIR}}/manta/{wildcards.sample} "\
		"--callRegions data/hs37d5_chr.bed.gz --exome;"\
		"python ${{TMP}}/manta/{wildcards.sample}/runWorkflow.py -j {params.threads} -g 16; "\
		"cp ${{TMPDIR}}/manta/{wildcards.sample}/results/variants/somaticSV.vcf.gz {output}.gz; "\
		"gunzip -d {output}.gz; "\


rule postprocess_manta_nocontrol:
	params:
		threads="1",
		runtime="3:45",
		memory="16000"
	input:
		VCF = config["working_dir"]+"/out_nocontrol/SV/manta/{sample}/{sample}_SV_unfiltered.vcf",
		BAM = config["working_dir"]+"/BAM/{sample}.bam",
		CNAs=config["working_dir"]+"/out_nocontrol/CNA/FREEC/{sample}/{sample}.bam_CNVs"
	output:
		config["working_dir"]+"/out_nocontrol/SV/manta/{sample}/{sample}_SV_filtered.vcf"
	conda:
		"../envs/WGS.yaml"
	shell:
		"python scripts/postprocess_mantaVCF.py -i {input.VCF} -o {output} "\
		"--minPR {config[manta_minPR]} --minSR {config[manta_minSR]} --minLen {config[manta_minLen_nocontrol]} "\
		"--pon {config[manta_pon]} --bam {input.BAM} --cnv {input.CNAs} --filterSmallInsertions {manta_filterSmallInsertions_nocontrol}"

rule postprocess_manta_control:
	params:
		threads="1",
		runtime="7:45",
		memory="16000"
	input:
		VCF = config["working_dir"]+"/out_control/SV/manta/{sample}/{sample}_SV_somatic.vcf",
		BAM = config["working_dir"]+"/BAM/{sample}.bam"
	output:
		config["working_dir"]+"/out_control/SV/manta/{sample}/{sample}_SV_filtered.vcf"
	conda:
		"../envs/WGS.yaml"
	shell:
		"python scripts/postprocess_mantaVCF.py -i {input.VCF} -o {output} "\
		"--minPR {config[manta_minPR]} --minSR {config[manta_minSR]} --minLen {config[manta_minLen_control]} "\
		"--bam {input.BAM} --tumorindex 1 --pon {config[manta_pon]} --filterSmallInsertions {manta_filterSmallInsertions_control}"

rule create_vizualization_datafiles_manta:
	params:
		threads="1",
		runtime="20",
		memory="8000"
	input:
		config["working_dir"]+"/{out}/SV/manta/{sample}/{sample}_SV_filtered.vcf"
	output:
		circos = config["working_dir"]+"/{out}/SV/manta/{sample}/{sample}_SV_circos.txt",
		chrplot = config["working_dir"]+"/{out}/SV/manta/{sample}/{sample}_SV_chrplots.tsv",
	conda:
		"../envs/WGS.yaml"
	shell:
		"python scripts/manta2circos.py -i {input} -o {output.circos};"\
		"python scripts/manta2chrplots.py -i {input} -o {output.chrplot};"
