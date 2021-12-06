## Create a count table from crispresso output that looks like:
## MappingSequence  A   B   C   D  ... other bin names

import os
import pandas as pd
import numpy as np
from scripts.make_count_table import *

def make_flat_table(samplesheet, outfile):
    allele_tbls = []
    for idx, row in samplesheet.iterrows():
        file = row['variantCountFile']

        if (os.path.exists(file)):
            allele_tbl = pd.read_csv(file, sep='\t')
            allele_tbl['SampleID'] = row['SampleID']
            allele_tbls.append(allele_tbl)

    flat = pd.concat(allele_tbls, axis='index', ignore_index=True)
    flat.to_csv(outfile, sep='\t', index=False, compression='gzip')


rule aggregate_variant_counts:
    input:
        lambda wildcards:
            samplesheet.at[wildcards.SampleID,'CRISPRessoDir']
    output:
        counts='results/variantCounts/{SampleID}.variantCounts.txt'
    run:
        aggregate_variant_counts(samplesheet, wildcards.SampleID, output.counts, config['variant_info'])


rule make_count_table_per_PCRrep:
    input:
        lambda wildcards:
            samplesheet.loc[samplesheet['ExperimentIDPCRRep']==wildcards.ExperimentIDPCRRep]['variantCountFile']
    output:
        counts='results/byPCRRep/{ExperimentIDPCRRep}.bin_counts.txt',
        freq='results/byPCRRep/{ExperimentIDPCRRep}.bin_freq.txt'
    run:
        make_count_table(samplesheet, 'ExperimentIDPCRRep', wildcards.ExperimentIDPCRRep, get_bin_list(), output.counts, output.freq)


rule make_count_table_per_experimentalRep:
    input:
        lambda wildcards:
            samplesheet.loc[samplesheet['ExperimentIDReplicates']==wildcards.ExperimentIDReplicates]['variantCountFile']
    output:
        counts='results/byExperimentRep/{ExperimentIDReplicates}.bin_counts.txt',
        freq='results/byExperimentRep/{ExperimentIDReplicates}.bin_freq.txt'
    run:
        make_count_table(samplesheet, 'ExperimentIDReplicates', wildcards.ExperimentIDReplicates, get_bin_list(), output.counts, output.freq, variantInfo=config['variant_info'])


rule write_pcr_replicate_correlation:
    input:
        variantCounts="results/summary/VariantCounts.DesiredVariants.flat.tsv",
        samplesheet="SampleList.snakemake.tsv"
    output:
        corfile="results/summary/PCRReplicateCorrelations.tsv",
        cvfile="results/summary/VariationVsAlleleFrequency.tsv",
        lowcorfile="results/summary/PCRReplicateCorrelations.LowQualSamples.tsv"
    params:
        codedir=config['codedir']
    shell:
        "Rscript {params.codedir}/workflow/scripts/GetPCRReplicateCorrelation.R \
          --variantCounts {input.variantCounts} \
          --samplesheet {input.samplesheet} \
          --correlationFile {output.corfile} \
          --lowCorSamplesFile {output.lowcorfile} \
          --cvFile {output.cvfile}"


rule make_count_table_per_experimentalRep_withCorFilter:
    input:
        samplesToExclude='results/summary/PCRReplicateCorrelations.LowQualSamples.tsv',
        samples = lambda wildcards:
            samplesheet.loc[samplesheet['ExperimentIDReplicates']==wildcards.ExperimentIDReplicates]['variantCountFile']
    output:
        counts='results/byExperimentRepCorFilter/{ExperimentIDReplicates}.bin_counts.txt',
        freq='results/byExperimentRepCorFilter/{ExperimentIDReplicates}.bin_freq.txt'
    run:
        make_count_table(samplesheet, 'ExperimentIDReplicates', wildcards.ExperimentIDReplicates, get_bin_list(), output.counts, output.freq, input.samplesToExclude)


rule make_flat_count_table:
    input:
        lambda wildcards: samplesheet['variantCountFile']
    output:
        flat='results/summary/VariantCounts.flat.tsv.gz'
    run:
        make_flat_table(samplesheet, output.flat) 


rule make_variant_matrix:
    input:
        variantCounts='results/summary/VariantCounts.flat.tsv.gz'
    output:
        matrix="results/summary/VariantCounts.matrix.tsv.gz"
    run:
        countsFlat = pd.read_csv(input.variantCounts, sep='\t')
        countsMatrix = countsFlat.pivot(
            index=[
                'MappingSequence',
                'AmpliconID',
                'VariantID',
                'RefAllele'],
            columns='SampleID',
            values='#Reads').reset_index()
        countsMatrix = countsMatrix.rename_axis(None, axis=1).reset_index(drop=True).fillna(0)  
        # the index needs to be reset because pivot names the index "SampleID"
        countsMatrix.to_csv(output.matrix, sep='\t', index=False)
