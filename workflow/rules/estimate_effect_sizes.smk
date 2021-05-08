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


# use the alleles from above to write a design file for feeding into mle
rule write_design_file:
	input: 
		counts='results/byExperimentRep/{ExperimentIDReplicates}.bin_counts.txt'
	output:
		design='results/byExperimentRep/{ExperimentIDReplicates}.design.txt'
	run:
		write_design_file(input.counts, output.design)


def get_sortparams_file(wildcards):
	currSamples = samplesheet.loc[(samplesheet['ExperimentIDReplicates'] == wildcards.ExperimentIDReplicates) & (samplesheet['Bin'].isin(binList))]
	Batch = currSamples['Batch'].unique()
	SampleNumber = currSamples['SampleNumber'].unique()
	if (len(Batch) != 1) or (len(SampleNumber) != 1):
		print(currSamples['SampleID'])
		raise ValueError("Found more than one possible sort params file path. Correct the samplesheet and rerun.")
	return os.path.join(config['sortparamsdir'], Batch[0] + "_" + SampleNumber[0] + ".csv")


# run mle to calculate the effect sizes
rule calculate_allelic_effect_sizes:
	input:
		counts='results/byExperimentRep/{ExperimentIDReplicates}.bin_counts.txt',
		design='results/byExperimentRep/{ExperimentIDReplicates}.design.txt',
		sortparams=get_sortparams_file 
	output:
		'results/byExperimentRep/{ExperimentIDReplicates}.raw_effects.txt'
	log:
		'results/byExperimentRep/{ExperimentIDReplicates}.mle_log.txt'
	shell:
		"""
		Rscript variant-flowfish/workflow/scripts/get_allele_effect_sizes.R \
			 --designDocLocation {input.design} \
			 --countsLocation {input.counts} \
			 --sortParamsloc {input.sortparams} \
			 --outputmle {output} --log {log}
		"""
