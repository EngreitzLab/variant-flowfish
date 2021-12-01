## Estimate effect sizes of each variant by calling the MLE on bin data


def write_design_file(counts_file, design):
	counts = pd.read_table(counts_file, index_col=0)
	ref_seq = counts[counts.columns[0]].idxmax() # note, this may not always work
	df = pd.DataFrame()
	df['MappingSequence'] = counts.index
	df['OligoID'] = df.MappingSequence
	df['OffTargetScore'] = 200
	df['chr'] = 'NA'
	df['start'] = 0
	df['end'] = 0
	df['name'] = '*'
	df['score'] = '*'
	df['strand'] = '*'
	df['GuideSequence'] = '*'
	df['GuideSequenceMinusG'] = '*'
	df['target'] = '*'
	df['subpool'] = '*'
	design_doc = df[['chr', 'start', 'end', 'name', 'score', 'strand', 'GuideSequence', 'GuideSequenceMinusG', 'MappingSequence', 'OffTargetScore', 'target', 'subpool', 'OligoID']]
	design_doc.loc[design_doc['MappingSequence'] == ref_seq, 'target'] = 'negative_control'
	design_doc.to_csv(design, sep='\t')


def get_sortparams_file(wildcards):
	if wildcards.directory == "byExperimentRep" or wildcards.directory == "byExperimentRepCorFilter":
		currSamples = samplesheet.loc[(samplesheet['ExperimentIDReplicates'] == wildcards.ExperimentID) & (samplesheet['sortParamsFile'] != "")]
	else:
		currSamples = samplesheet.loc[(samplesheet['ExperimentIDPCRRep'] == wildcards.ExperimentID) & (samplesheet['sortParamsFile'] != "")]
	sortParams = currSamples['sortParamsFile'].unique()
	if (len(sortParams) > 1):
		print(sortParams)
		print(wildcards)
		raise ValueError("Found more than one possible sort params file path. Correct the samplesheet and rerun.")
	elif (len(sortParams) == 0):
		print(sortParams)
		print(wildcards)
		raise ValueError("Found no possible sort params file path. Correct the samplesheet and rerun.")
	else:
		return sortParams[0]


# run mle to calculate the effect sizes
rule calculate_allelic_effect_sizes:
	input:
		counts='results/{directory}/{ExperimentID}.bin_counts.txt',
		sortparams=get_sortparams_file 
	output:
		'results/{directory}/{ExperimentID}.raw_effects.txt'
	log:
		'results/{directory}/{ExperimentID}.mle_log.txt'
	shell:
		"""
		Rscript variant-flowfish/workflow/scripts/get_allele_effect_sizes.R \
			 --countsLocation {input.counts} \
			 --sortParamsloc {input.sortparams} \
			 --outputmle {output} --log {log}
		"""


# Normalize allele effect sizes to reference allele specified in the variant info table
rule normalize_allelic_effect_sizes:
	input:
		'results/{directory}/{ExperimentID}.raw_effects.txt'
	output:
		'results/{directory}/{ExperimentID}.effects_vs_ref.txt'
	params:
		codedir=config['codedir'],
		variantInfo=config['variant_info']
	shell:
		"""
		python {params.codedir}/workflow/scripts/normalize_allele_effects.py -i {input} -o {output} -v {params.variantInfo}
		"""

rule plot_allelic_effect_sizes:
	input:
		'results/{replicateDirectory}/{ExperimentIDReplicates}.effects_vs_ref.txt'
	output:
		'results/{replicateDirectory}/{ExperimentIDReplicates}.effects_vs_ref.pdf'
	params:
		codedir=config['codedir']
	shell:
		"""
		Rscript {params.codedir}/workflow/scripts/PlotMleVariantEffects.R --mleEffects {input} --outfile {output} 
		"""


###################################################################################
## Run again, but ignore input bin counts

# run mle to calculate the effect sizes
rule calculate_allelic_effect_sizes_ignoreInputBin:
	input:
		counts='results/{directory}/{ExperimentID}.bin_counts.topN.txt',
		sortparams=get_sortparams_file 
	output:
		'results/{directory}/{ExperimentID}.raw_effects_ignoreInputBin.txt'
	log:
		'results/{directory}/{ExperimentID}.mle_log_ignoreInputBin.txt'
	shell:
		"""
		Rscript variant-flowfish/workflow/scripts/get_allele_effect_sizes.R \
			 --countsLocation {input.counts} \
			 --sortParamsloc {input.sortparams} \
			 --outputmle {output} --log {log} \
			 --ignoreInputBinCounts TRUE
		"""

# Normalize allele effect sizes to reference allele specified in the variant info table
rule normalize_allelic_effect_sizes_ignoreInputBin:
	input:
		'results/{directory}/{ExperimentID}.raw_effects_ignoreInputBin.txt'
	output:
		'results/{directory}/{ExperimentID}.effects_vs_ref_ignoreInputBin.txt'
	params:
		codedir=config['codedir'],
		variantInfo=config['variant_info']
	shell:
		"""
		python {params.codedir}/workflow/scripts/normalize_allele_effects.py -i {input} -o {output} -v {params.variantInfo}
		"""

rule plot_allelic_effect_sizes_ignoreInputBin:
	input:
		'results/{replicateDirectory}/{ExperimentIDReplicates}.effects_vs_ref_ignoreInputBin.txt'
	output:
		'results/{replicateDirectory}/{ExperimentIDReplicates}.effects_vs_ref_ignoreInputBin.pdf'
	params:
		codedir=config['codedir']
	shell:
		"""
		Rscript {params.codedir}/workflow/scripts/PlotMleVariantEffects.R --mleEffects {input} --outfile {output} 
		"""
