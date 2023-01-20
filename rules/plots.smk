rule make_linearplots:
	params:
		threads="1",
		runtime="4:00",
		memory="16000",
		sex = lambda wcs: sample_sex[wcs.sample]
	input:
		CNA = config["output_dir"]+"/out/CNA/FREEC/{sample}/{sample}_CNA_chrplots.tsv",
		SV = config["output_dir"]+"/out/SV/manta/{sample}/{sample}_SV_chrplots.tsv"
	output:
		expand(config["output_dir"]+"/out/plots/chrplots_{{format}}/{{sample}}/{{sample}}_chr{chr}.{{format}}",chr=chromosomes)
	shell:
		"python scripts/plot_chr_CNA-SV.py --cnv {input.CNA} --sv {input.SV} --sample {wildcards.sample} --sex {params.sex} -o {config[output_dir]}/out/plots/chrplots_{wildcards.format}/{wildcards.sample}/{wildcards.sample} --format {wildcards.format}"
