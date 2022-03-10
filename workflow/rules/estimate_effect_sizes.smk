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
		raw='results/{directory}/{ExperimentID}.raw_effects.txt',
		counts='results/{directory}/{ExperimentID}.rescaled_counts.txt'	
	log:
		'results/{directory}/{ExperimentID}.mle_log.txt'
	shell:
		"""
		Rscript variant-flowfish/workflow/scripts/get_allele_effect_sizes.R \
			 --countsLocation {input.counts} \
			 --sortParamsloc {input.sortparams} \
			 --outputmle {output.raw} --log {log} --rescaledCounts {output.counts}
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


def aggregate_allelic_effect_tables(samplesheet, outfile, wildcards):
    allele_tbls = []
    for expt in samplesheet[wildcards.ExperimentID].drop_duplicates():
        file = 'results/{dir}/{e}.effects_vs_ref_ignoreInputBin.txt'.format(dir=wildcards.replicateDirectory, e=expt)

        if (os.path.exists(file)):
            allele_tbl = pd.read_csv(file, sep='\t')
            allele_tbl[wildcards.ExperimentID] = expt
            id_split = allele_tbl[wildcards.ExperimentID].str.split('-', expand=False)
            if wildcards.ExperimentID == "ExperimentIDReplicates": # byExperimentRep
                allele_tbl['BioRep'] = id_split.str[-3]
                allele_tbl['Spike-In'] = id_split.str[-2]
                allele_tbl['FFRep'] = id_split.str[-1]
            else: #byPCRRep
                assert wildcards.ExperimentID == "ExperimentIDPCRRep"
                allele_tbl['BioRep'] = id_split.str[-4]
                allele_tbl['Spike-In'] = id_split.str[-3]
                allele_tbl['FFRep'] = id_split.str[-2]
                allele_tbl['PCRRep'] = id_split.str[-1]
            allele_tbl['Location'] = allele_tbl['VariantID'].str.split(':').str[:3].str.join(':')
            allele_tbl['Variant'] = allele_tbl['VariantID'].str.split(':').str[-1:].str.join(':')  
            allele_tbls.append(allele_tbl)

    flat = pd.concat(allele_tbls, axis='index', ignore_index=True)
    flat.to_csv(outfile, sep='\t', index=False, compression='gzip')


rule aggregate_allelic_effects:
    input:
        lambda wildcards: 
        	['results/{dir}/{e}.effects_vs_ref_ignoreInputBin.txt'.format(dir=wildcards.replicateDirectory, e=e) for e in samplesheet.loc[samplesheet['Bin'].isin(binList),wildcards.ExperimentID].drop_duplicates()]
    output:
        flat='results/summary/AllelicEffects.{replicateDirectory}.{ExperimentID}.flat.tsv.gz'
    run:
        aggregate_allelic_effect_tables(samplesheet, output.flat, wildcards) 


###################################################################################
## Run again, but ignore input bin counts

# run mle to calculate the effect sizes
rule calculate_allelic_effect_sizes_ignoreInputBin:
	input:
		counts='results/{directory}/{ExperimentID}.bin_counts.txt',
		sortparams=get_sortparams_file 
	output:
		raw='results/{directory}/{ExperimentID}.raw_effects_ignoreInputBin.txt',
		counts='results/{directory}/{ExperimentID}.rescaled_counts_ignoreInputBin.txt'	
	log:
		'results/{directory}/{ExperimentID}.mle_log_ignoreInputBin.txt'
	shell:
		"""
		Rscript variant-flowfish/workflow/scripts/get_allele_effect_sizes.R \
			 --countsLocation {input.counts} \
			 --sortParamsloc {input.sortparams} \
			 --outputmle {output.raw} --log {log} --rescaledCounts {output.counts} \
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
