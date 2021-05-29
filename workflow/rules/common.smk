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

	for read in ["1","2"]:
		colName = 'fastqR' + read
		if not colName in samplesheet.columns:
			samplesheet[colName] = ""
			for i in samplesheet.index:
				currSample = samplesheet.at[i,'SampleID']
				file = glob.glob("{}_*_R{}_*fastq.gz".format(os.path.join(fastqdir, currSample), read))
				if len(file) > 1:
					raise ValueError("Found more than one FASTQ file for sample :" + currSample)
				if len(file) == 1:
					samplesheet.at[i,colName] = file[0]
	
	return samplesheet


def add_experiment_names(samplesheet):
	if ('ExperimentIDReplicates' in samplesheet.columns) or ('ExperimentID' in samplesheet.columns) or ('ExperimentIDPCRRep' in samplesheet.columns):
		print("Warning: ExperimentID columns found and will be overwritten in the sample sheet")

	## Experiments at the level of PCR replicates 
	s = samplesheet[keyCols + ['VFFSpikeIn'] + repCols + ['PCRRep']].drop_duplicates()
	s['ExperimentIDPCRRep'] = ['-'.join([str(v) for v in list(row.values)]) for index,row in s.iterrows()]
	samplesheet = samplesheet.merge(s)

	## Experiments at the level of specified replicate columns, before spike correction
	s = samplesheet[keyCols + ['VFFSpikeIn'] + repCols].drop_duplicates()
	s['ExperimentIDReplicates'] = ['-'.join([str(v) for v in list(row.values)]) for index,row in s.iterrows()]
	samplesheet = samplesheet.merge(s)

	## Experiments at the level of experiments (combined across replicates)
	s = samplesheet[keyCols].drop_duplicates()
	s['ExperimentID'] = ['-'.join([str(v) for v in list(row.values)]) for index,row in s.iterrows()]
	samplesheet = samplesheet.merge(s)

	return(samplesheet)	


def add_outputs(samplesheet):
	samplesheet['CRISPRessoDir'] = ['results/crispresso/CRISPResso_on_{SampleID}/'.format(SampleID=row['SampleID']) for index, row in samplesheet.iterrows()]
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
	samplesheet = add_experiment_names(samplesheet)
	samplesheet = add_outputs(samplesheet)
	samplesheet.index = samplesheet['SampleID']

	return samplesheet


def get_bin_list():
	binList = samplesheet['Bin'].drop_duplicates()
	binList = binList[(binList != "All") & (binList != "Neg") & (binList.notnull())]
	binList = [str(b) for b in list(binList)]
	print("Processing unique bins: " + ' '.join(binList))
	return(binList)



# global variables
requiredCols = ['SampleID','AmpliconID','Bin','PCRRep','VFFSpikeIn']
ampliconRequiredCols = ['AmpliconID','AmpliconSeq','GuideSpacer']  ## To do:  Allow specifying crispresso quantification window for different amplicons
keyCols = config['experiment_keycols'].split(',')
repCols = config['replicate_keycols'].split(',')
codedir = config['codedir']
sortparamsdir = config['sortparamsdir']
fastqdir = config['fastqdir']
#n_reps = 2 # config['n_reps']

samplesheet = load_sample_sheet(config['sample_sheet'], config['amplicon_info'])
samplesheet.to_csv("SampleList.snakemake.tsv", index=False, header=True, sep='\t')
binList = get_bin_list()



#######################################################################################
####### helpers ###########

def all_input(wildcards):

	wanted_input = []

	## CRISPResso output:
	wanted_input.extend(list(samplesheet['CRISPRessoDir'].unique()))

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

	## Output files for PCR replicates (before merging spike-in data)
	wanted_input.extend(list(samplesheet['ExperimentIDPCRRep_BinCounts'].unique()))
	wanted_input.append("results/summary/VariantCounts.flat.tsv.gz")

	## Output files for PCR replicates (after merging spike-in data) (?)
	wanted_input.extend([])

	## Output files for replicate experiments (before merging spike-in data)
	if len(repCols) > 0:
		wanted_input.extend(list(samplesheet['ExperimentIDReplicates_BinCounts'].unique()))
		wanted_input.extend([
			'results/byExperimentRep/{}.raw_effects.txt'.format(e) for e in samplesheet.loc[samplesheet['Bin'].isin(binList)]['ExperimentIDReplicates'].unique()
		])


	## Output files for replicate experiments (after merging spike-in data)
	wanted_input.extend([])

	## Output files for experiments (with replicates merged, before merging spike-in data)
	wanted_input.extend([])

	## Output files for experiments (with replicates merged, after merging spike-in data)
	wanted_input.extend([])

	return wanted_input


