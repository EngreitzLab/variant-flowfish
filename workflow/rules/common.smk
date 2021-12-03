from snakemake.utils import validate
import pandas as pd
import os
import glob

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
# singularity: "docker://continuumio/miniconda3"



###########################################################################################
##### load config and sample sheets #####

# configfile: "config/config.yaml"   ## Read from command line instead
# validate(config, schema="../schemas/config.schema.yaml")


def find_fastq_files(samplesheet, fastqdir):
	## Adds columns 'fastqR1' and 'fastqR2' to the sample sheet, only if they do not already exist

	if single_end:
		reads = ["1"]
	else:
		reads = ["1","2"]

	for read in reads:
		colName = 'fastqR' + read
		if not colName in samplesheet.columns:
			samplesheet[colName] = ""
			for i in samplesheet.index:
				currSample = samplesheet.at[i,'SampleID']
				## Try looking for the bcl2fastq output file format
				file = glob.glob("{}_*_R{}_*fastq.gz".format(os.path.join(fastqdir, currSample), read))
				if len(file) == 0: ## Try looking for the barcode-splitter output format
					file = glob.glob("{}-read-{}.fastq.gz".format(os.path.join(fastqdir, currSample), read))
				if len(file) > 1:
					raise ValueError("Found more than one FASTQ file for sample :" + currSample)
				elif len(file) == 0:
					print("Warning: Could not find FASTQ file for read " + read + " and sample: " + currSample)
				elif len(file) == 1:
					samplesheet.at[i,colName] = file[0]

	
	return samplesheet


def find_sort_params_files(samplesheet):
	## If the user provided 'sortParamsFile' in the samplesheet, use that. Otherwise, look in the sortParams directory
	if not 'sortParamsFile' in samplesheet.columns:
		samplesheet['sortParamsFile'] = [os.path.join(config['sortparamsdir'], str(row['Batch']) + "_" + str(row['SampleNumber']) + ".csv") for idx, row in samplesheet.iterrows()]
		samplesheet.loc[samplesheet['SampleNumber'].isnull(),'sortParamsFile'] = ""
	return samplesheet


def add_experiment_names(samplesheet):
	if ('ExperimentIDReplicates' in samplesheet.columns) or ('ExperimentID' in samplesheet.columns) or ('ExperimentIDPCRRep' in samplesheet.columns):
		print("Warning: ExperimentID columns found and will be overwritten in the sample sheet")

	## Experiments at the level of PCR replicates 
	if genotyping_only:
		cols_ExperimentIDPCRRep = keyCols + repCols + ['PCRRep']
		cols_ExperimentIDReplicates = keyCols + repCols
	else:
		cols_ExperimentIDPCRRep = keyCols + ['VFFSpikeIn'] + repCols + ['PCRRep']
		cols_ExperimentIDReplicates = keyCols + ['VFFSpikeIn'] + repCols

	s = samplesheet[cols_ExperimentIDPCRRep].drop_duplicates()
	s['ExperimentIDPCRRep'] = ['-'.join([str(v) for v in list(row.values)]) for index,row in s.iterrows()]
	samplesheet = samplesheet.merge(s)

	## Experiments at the level of specified replicate columns, before spike correction
	s = samplesheet[cols_ExperimentIDReplicates].drop_duplicates()
	s['ExperimentIDReplicates'] = ['-'.join([str(v) for v in list(row.values)]) for index,row in s.iterrows()]
	samplesheet = samplesheet.merge(s)

	## Experiments at the level of experiments (combined across replicates)
	s = samplesheet[keyCols].drop_duplicates()
	s['ExperimentID'] = ['-'.join([str(v) for v in list(row.values)]) for index,row in s.iterrows()]
	samplesheet = samplesheet.merge(s)

	return(samplesheet)	


def add_outputs(samplesheet):
	samplesheet['CRISPRessoDir'] = ['results/crispresso/CRISPResso_on_{SampleID}/'.format(SampleID=row['SampleID']) for index, row in samplesheet.iterrows()]
	if not genotyping_only:
		samplesheet['ExperimentIDPCRRep_BinCounts'] = ['results/byPCRRep/{}.bin_counts.txt'.format(e) for e in samplesheet['ExperimentIDPCRRep']]
		samplesheet['ExperimentIDReplicates_BinCounts'] = ['results/byExperimentRep/{}.bin_counts.txt'.format(e) for e in samplesheet['ExperimentIDReplicates']]
	return samplesheet


def validate_sample_sheet(samplesheet):
	print("Validating the Sample Sheet ...\n")

	if not samplesheet['SampleID'].is_unique:
		raise ValueError("SampleID column in samplesheet must not contain duplicates.")

	for col in requiredCols:
		if not col in samplesheet.columns:
			raise ValueError("Missing required column in sample sheet: " + col)

	for col in keyCols:
		if not col in samplesheet.columns:
			raise ValueError("Missing column in sample sheet that is provided in experiment_keycols in the config file: " + col)

	for col in repCols:
		if not col in samplesheet.columns:
			raise ValueError("Missing column in sample sheet that is provided in replicate_keycols in the config file: " + col)

	print("Found all Experiment Key columns. Generating comparisons for the following experiments:")
	print('\t'.join(keyCols))
	for index, row in samplesheet[keyCols].drop_duplicates().iterrows():
		print('\t'.join([str(v) for v in row.values]))


def load_sample_sheet(samplesheetFile, ampliconInfoFile, idcol='AmpliconID'):
	samplesheet = pd.read_table(samplesheetFile, dtype=str)
	validate_sample_sheet(samplesheet)
	samplesheet.index = samplesheet['SampleID']  ## Requires that SampleID is unique

	if not set(ampliconRequiredCols).issubset(samplesheet.columns):
		amplicons = pd.read_table(ampliconInfoFile)
		if not set(samplesheet[idcol]).issubset(amplicons[idcol]):
			raise ValueError("Some AmpliconIDs in the sample sheet are not specified in the amplicon info table.")
		if not set(ampliconRequiredCols).issubset(amplicons.columns):
			raise ValueError("Amplicon info file must contain AmpliconID AmpliconSeq GuideSpacer")
		samplesheet = samplesheet.merge(amplicons[ampliconRequiredCols])
		if not set(ampliconRequiredCols).issubset(samplesheet.columns):
			raise ValueError("Failed to merge samplesheet and amplicon info file.")		

	samplesheet = find_fastq_files(samplesheet, fastqdir)
	if not genotyping_only:
		samplesheet = find_sort_params_files(samplesheet)
	samplesheet = add_experiment_names(samplesheet)
	samplesheet = add_outputs(samplesheet)
	samplesheet.index = samplesheet['SampleID']

	return samplesheet


def get_bin_list():
	binList = samplesheet['Bin'].drop_duplicates()
	if  "All" not in binList.tolist():
		print("\nWARNING: Did not find any entries with Bin == 'All' (unsorted edited cells input into FlowFISH). Was this intended, or was 'All' mispelled?\n\n")
	if "Neg" not in binList.tolist():
		print("\nWARNING: Did not find any entries with Bin == 'Neg' (unedited cells used to assess sequencing error rate). Was this intended, or was 'Neg' mispelled?\n\n")		
	binList = binList[(binList != "All") & (binList != "Neg") & (binList.notnull())]
	binList = [str(b) for b in list(binList)]
	print("Processing unique bins: " + ' '.join(binList))
	return(binList)

def load_variant_Table(variant_table, requiredCols):
	variants = pd.read_table(variant_table, dtype=str)
	if not set(requiredCols).issubset(variants.columns):
		raise ValueError("Variant table is missing required cols: ", set(requiredCols).difference(variants.columns))
	return(variants)


# global variables
genotyping_only = ('genotyping_only' in config) and (bool(config['genotyping_only']))
if genotyping_only:
	requiredCols = ['SampleID','AmpliconID','Bin','PCRRep']
else:
	requiredCols = ['SampleID','AmpliconID','Bin','PCRRep','VFFSpikeIn']

single_end = ('single_end' in config) and (bool(config['single_end']))

ampliconRequiredCols = ['AmpliconID','AmpliconSeq','GuideSpacer']  ## To do:  Allow specifying crispresso quantification window for different amplicons
variantRequiredCols = ['AmpliconID','VariantID','MappingSequence','RefAllele']
keyCols = config['experiment_keycols'].split(',')
repCols = config['replicate_keycols'].split(',')
codedir = config['codedir']
fastqdir = config['fastqdir']
sortparamsdir = config['sortparamsdir'] if not genotyping_only else None
#n_reps = 2 # config['n_reps']

samplesheet = load_sample_sheet(config['sample_sheet'], config['amplicon_info'])
samplesheet.to_csv("SampleList.snakemake.tsv", index=False, header=True, sep='\t')
binList = get_bin_list()

if 'variant_table' in config:
	variants = load_variant_table(config['variant_table'], variantRequiredCols)
else:
	variants = None


#######################################################################################
####### helpers ###########

def all_input(wildcards):

	wanted_input = []

	## CRISPResso output:
	wanted_input.extend(list(samplesheet['CRISPRessoDir'].unique()))
	wanted_input.append("results/summary/VariantCounts.flat.tsv.gz")
	wanted_input.append("results/summary/VariantCounts.DesiredVariants.flat.tsv")
	wanted_input.append("results/summary/VariantCounts.DesiredVariants.matrix.tsv")

	## Genotyping plots:
	wanted_input.append("results/summary/DesiredVariants.RData")

	#wanted_input.extend(
	#	['crispresso/CRISPResso_on_{SampleID}/{AmpliconID}.Alleles_frequency_table_around_sgRNA_{GuideSpacer}.txt'.format(
	#		SampleID=row['SampleID'], 
	#		AmpliconID=row['AmpliconID'], 
	#		GuideSpacer=row['GuideSpacer']) 
	#	for index, row in samplesheet.iterrows()]
	#)

	## Bowtie2 alignments:
	wanted_input.extend(
		['results/aligned/{s}/{s}.bam'.format(s=s) for s in samplesheet['SampleID'].unique()]
	)
	wanted_input.append("results/summary/alignment.counts.tsv")

	## At what point do we merge in the spike-in data? 

	if not genotyping_only:
		## Output files for PCR replicates (before merging spike-in data)
		wanted_input.extend(list(samplesheet['ExperimentIDPCRRep_BinCounts'].unique()))
		#wanted_input.extend([
		#	'results/byPcrRep/{}.effects_vs_ref.pdf'.format(e) for e in samplesheet.loc[samplesheet['Bin'].isin(binList)]['ExperimentIDPCRRep'].unique()
		#])

		## Output files for PCR replicates (after merging spike-in data) (?)
		wanted_input.extend([])

		## Output files for replicate experiments (before merging spike-in data)
		wanted_input.extend(list(samplesheet['ExperimentIDReplicates_BinCounts'].unique()))
		wanted_input.extend([
			'results/byExperimentRep/{}.effects_vs_ref.pdf'.format(e) for e in samplesheet.loc[samplesheet['Bin'].isin(binList)]['ExperimentIDReplicates'].unique()
		])

		wanted_input.extend([
			'results/byExperimentRepCorFilter/{}.effects_vs_ref.pdf'.format(e) for e in samplesheet.loc[samplesheet['Bin'].isin(binList)]['ExperimentIDReplicates'].unique()
		])
		wanted_input.extend([
			'results/byExperimentRepCorFilter/{}.effects_vs_ref_ignoreInputBin.pdf'.format(e) for e in samplesheet.loc[samplesheet['Bin'].isin(binList)]['ExperimentIDReplicates'].unique()
		])

		## Output files for replicate experiments (after merging spike-in data)
		wanted_input.extend([])

		## Output files for experiments (with replicates merged, before merging spike-in data)
		wanted_input.extend([])

		## Output files for experiments (with replicates merged, after merging spike-in data)
		wanted_input.extend([])

	return wanted_input


