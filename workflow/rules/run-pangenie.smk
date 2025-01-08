

rule pangenie_index:
	"""
	Create index for PanGenie.
	"""
	input:
		vcf = PANEL_MULTI,
		fasta = REFERENCE,
	output:
		directory("{results}/pangenie/index/")
	log:
		"{results}/pangenie/index/index.log"
	resources:
		mem_mb = 80000,
		walltime = "5:00:00"
	threads: 24
	params:
		out_prefix = "{results}/pangenie/index/index"
	benchmark:
		"{results}/pangenie/index/index_benchmark.txt"
	singularity:
		"workflow/container/eblerjana_eblerjana_pangenie-sampler.sif"
	shell:
		"""
		PanGenie-index -v {input.vcf} -r {input.fasta} -o {params.out_prefix} -t {threads}  &> {log}
		"""


rule pangenie_genotype_full:
	"""
	Run genotyping using the full panel (no sampling).
	"""
	input:
		reads = lambda wildcards: READS[wildcards.sample],
		index = directory("{results}/pangenie/index/")
	output:
		"{results}/pangenie/sampled-full/{sample}/{sample}-pangenie_multi_genotyping.vcf"
	log:
		"{results}/pangenie/sampled-full/{sample}/{sample}-pangenie_multi_genotyping.log"
	resources:
		mem_mb = 70000,
		walltime = "5:00:00"
	params:
		index = "{results}/pangenie/index/index",
		out_prefix = "{results}/pangenie/sampled-full/{sample}/{sample}-pangenie_multi"
	benchmark:
		"{results}/pangenie/sampled-full/{sample}/{sample}-pangenie_multi_benchmark.txt"
	threads:
		24
	singularity:
		"workflow/container/eblerjana_eblerjana_pangenie-sampler.sif"
	shell:
		"""
		PanGenie -f {params.index} -i <(gunzip -c {input.reads}) -o {params.out_prefix} -t {threads} -j {threads} -a 108 &> {log}
		"""


rule pangenie_genotype_sampling:
	"""
	Run genotyping with subsampling.
	"""
	input:
		reads = lambda wildcards: READS[wildcards.sample],
		index = directory("{results}/pangenie/index/")
	output:
		"{results}/pangenie/sampled-{size}/{sample}/{sample}-pangenie_multi_genotyping.vcf"
	log:
		"{results}/pangenie/sampled-{size}/{sample}/{sample}-pangenie_multi_genotyping.log"
	resources:
		mem_mb = 70000,
		walltime = "12:00:00"
	wildcard_constraints:
		size = "[0-9]+"
	params:
		index = "{results}/pangenie/index/index",
		out_prefix = "{results}/pangenie/sampled-{size}/{sample}/{sample}-pangenie_multi"
	benchmark:
		"{results}/pangenie/sampled-{size}/{sample}/{sample}-pangenie_multi_benchmark.txt"
	threads:
		24
	singularity:
		"workflow/container/eblerjana_eblerjana_pangenie-sampler.sif"
	shell:
		"""
		PanGenie -f {params.index} -i <(gunzip -c {input.reads}) -o {params.out_prefix} -t {threads} -j {threads} -d -x {wildcards.size} &> {log}
		"""


rule pangenie_sampling:
	"""
	Run just sampling.
	"""
	input:
		reads = lambda wildcards: READS[wildcards.sample],
		index = directory("{results}/pangenie/index/")
	output:
		"{results}/sampling/sampled-{size}/{sample}/{sample}-pangenie_multi_panel.vcf"
	log:
		"{results}/sampling/sampled-{size}/{sample}/{sample}-pangenie_multi_sampling.log"
	resources:
		mem_mb = 70000,
		walltime = "12:00:00"
	params:
		index = "{results}/pangenie/index/index",
		out_prefix = "{results}/sampling/sampled-{size}/{sample}/{sample}-pangenie_multi"
	benchmark:
		"{results}/sampling/sampled-{size}/{sample}/{sample}-pangenie_multi_benchmark.txt"
	threads:
		24
	singularity:
		"workflow/container/eblerjana_eblerjana_pangenie-sampler.sif"
	shell:
		"""
		PanGenie-sampling -f {params.index} -i <(gunzip -c {input.reads}) -o {params.out_prefix} -t {threads} -j {threads} -x {wildcards.size} &> {log}
		"""


rule convert_to_biallelic:
	"""
	Convert the bubble genotypes to variant genotypes.
	(Only possible if necessary annotations are present)
	"""
	input:
		genotypes = "{results}/pangenie/sampled-{size}/{sample}/{sample}-pangenie_multi_genotyping.vcf",
		panel_bi = PANEL_BI
	output:
		"{results}/pangenie/sampled-{size}/{sample}/{sample}-pangenie_bi_genotyping.vcf.gz"
	resources:
		mem_mb = 20000,
		walltime = "2:00:00"
	shell:
		"""
		cat {input.genotypes} | python3 workflow/scripts/convert-to-biallelic.py {input.panel_bi} | bgzip > {output} 
		"""
