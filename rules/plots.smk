rule make_chrplots:
	params:
		threads="1",
		runtime="4:00",
		memory="16000",
		sex = lambda wcs: df_samplesheet.loc[wcs.sample,"sex"]
    conda:
		"../envs/WGS.yaml"
	input:
		CNA = config["working_dir"]+"/{out}/CNA/FREEC/{sample}/{sample}_CNA_chrplots.tsv",
		SV = config["working_dir"]+"/{out}/SV/manta/{sample}/{sample}_SV_chrplots.tsv"
	output:
		expand(config["working_dir"]+"/{{out}}/plots/chrplots_{{format}}/{{sample}}/{{sample}}_chr{chr}.{{format}}",chr=chromosomes)
	shell:
		"python scripts/plot_chr_CNA-SV.py --cna {input.CNA} --sv {input.SV} --sample {wildcards.sample} --sex {params.sex} -o {config[working_dir]}/{wildcards.out}/plots/chrplots_{wildcards.format}/{wildcards.sample}/{wildcards.sample} --format {wildcards.format}"


rule make_circosplots:
	params:
		threads="1",
		runtime="4:00",
		memory="16000",
		sex = lambda wcs: df_samplesheet.loc[wcs.sample,"sex"]
    conda:
		"../envs/WGS.yaml"
	input:
		CNA = config["working_dir"]+"/{out}/CNA/FREEC/{sample}/{sample}_CNA_chrplots.tsv",
		SV = config["working_dir"]+"/{out}/SV/manta/{sample}/{sample}_SV_circos.txt"
	output:
		config["working_dir"]+"/{out}/plots/circos/circos_{sample}.svg"
	shell:
		"python scripts/plot_circos.py --cna {input.CNA} --sv {input.SV} --sex {params.sex} --data data -o {output}"
