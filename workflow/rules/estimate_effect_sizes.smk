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
		counts='{path}.bin_counts.topN.txt',
		sortparams=get_sortparams_file 
	output:
		'{path}.raw_effects.txt'
	log:
		'{path}.mle_log.txt'
	shell:
		"""
		Rscript variant-flowfish/workflow/scripts/get_allele_effect_sizes.R \
			 --countsLocation {input.counts} \
			 --sortParamsloc {input.sortparams} \
			 --outputmle {output} --log {log}
		"""






# Rscript variant-flowfish/workflow/scripts/get_allele_effect_sizes.R \
# 			 --countsLocation '/oak/stanford/groups/engreitz/Projects/VariantEditing/FF/210426B001-VFF/results/byExperimentRep/IL2RA-Jurkat-IL2RA_peg10-0-1.bin_counts.txt' \
# 			 --sortParamsloc '/oak/stanford/groups/engreitz/Projects/VariantEditing/FF/210426B001-VFF/sortParams/tmp/210426B001_3.csv' \
# 			 --outputmle '/oak/stanford/groups/engreitz/Projects/VariantEditing/FF/210426B001-VFF/test/test.raw_effects.txt' --log '/oak/stanford/groups/engreitz/Projects/VariantEditing/FF/210426B001-VFF/test/test.mle_log.txt'