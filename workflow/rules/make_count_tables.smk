## Create a count table from crispresso output that looks like:
## MappingSequence  A   B   C   D  ... other bin names

import os
import pandas as pd
import numpy as np

def make_count_table(samplesheet, group_col, group_id, bins, outfile, outfile_frequencies, samplesToExclude=None, edit_regions=None):
    ## Function to make a count table at various layers of resolution (e.g., by experiment, or by replicate, or by PCR replicate)
    ## To do: Move the python code for these rules into separate python scripts so they can be run independently of the snakemake pipeline (at least, this makes it easier to test and debug the code)

    currSamples = samplesheet.loc[samplesheet[group_col]==group_id]

    if samplesToExclude is not None:
        exclude = pd.read_table(samplesToExclude, header=None, names=['SampleID'])
        currSamples = currSamples[~currSamples['SampleID'].isin(exclude['SampleID'].values)]

    allele_tbls = []

    for idx, row in currSamples.iterrows():
        file = "results/crispresso/CRISPResso_on_{SampleID}/Alleles_frequency_table.zip".format(
            SampleID=row['SampleID'])

        if (os.path.exists(file)):
            # import pdb; pdb.set_trace()
            allele_tbl = pd.read_csv(file, sep='\t')
            allele_tbl['#Reads'] = allele_tbl['#Reads'].astype(np.int32)
            # check only edit region for matching? 
            if edit_regions:
                sample_id = row['SampleID']
                amplicon = sample_id.split('-')[4] # this is not robust but will do for now
                amplicon = amplicon.lstrip('Amplicon')
                if amplicon not in edit_regions.keys():
                    print('Invalid amplicon name.')
                    exit(1)
                start, end = edit_regions[amplicon].split('-')
                start = int(start)
                end = int(end)
                allele_tbl['Aligned_Sequence_edit_region'] = allele_tbl['Aligned_Sequence'].apply(lambda x: x[start:end+1])
                allele_tbl['Reference_Sequence_edit_region'] = allele_tbl['Reference_Sequence'].apply(lambda x: x[start:end+1])
                ref_seq = allele_tbl[allele_tbl['Aligned_Sequence_edit_region'] == allele_tbl['Reference_Sequence_edit_region']]['Reference_Sequence'].mode().item()
            else:
                ref_seq = allele_tbl[allele_tbl['Aligned_Sequence'] == allele_tbl['Reference_Sequence']]['Reference_Sequence'].mode().item()
            allele_tbl = allele_tbl[allele_tbl['Reference_Sequence'] == ref_seq] # necessary?
            allele_tbl = allele_tbl[['Aligned_Sequence', '#Reads']]
            allele_tbl.columns = ['Aligned_Sequence', row['SampleID']]
            allele_tbls.append(allele_tbl)

    if len(allele_tbls) > 0:
        count_tbl = allele_tbls.pop()

        for tbl in allele_tbls:
            count_tbl = count_tbl.merge(tbl, on='Aligned_Sequence', how='outer')

        count_tbl = count_tbl.set_index('Aligned_Sequence')
        count_tbl = count_tbl.fillna(0)
        count_tbl = count_tbl.astype(int)

        ## Now, sum counts per bin
        bin_list = bins + list(set(currSamples['Bin'].unique())-set(bins))
        for uniqBin in bin_list:
            samples = currSamples.loc[currSamples['Bin'] == uniqBin]
            if len(samples) > 0:
                count_tbl[uniqBin] = count_tbl[samples['SampleID']].sum(axis=1).values

            else:
                count_tbl[uniqBin] = 0
        count_tbl = count_tbl[bin_list]

    else:
        count_tbl = pd.DataFrame({'Aligned_Sequence':[]})
        for uniqBin in bins:
            count_tbl[uniqBin] = []
        count_tbl = count_tbl.set_index('Aligned_Sequence')

    count_tbl.index.name = "MappingSequence"
    count_tbl.to_csv(outfile, sep='\t')

    freq_tbl = count_tbl.div(count_tbl.sum(axis=0), axis=1)
    freq_tbl.to_csv(outfile_frequencies, sep='\t', float_format='%.6f')



def make_flat_table(samplesheet, outfile):

    allele_tbls = []
    for idx, row in samplesheet.iterrows():
        file = "{CRISPRessoDir}/Alleles_frequency_table.zip".format(
            CRISPRessoDir=row['CRISPRessoDir'])

        if (os.path.exists(file)):
            allele_tbl = pd.read_csv(file, sep='\t')
            allele_tbl['SampleID'] = row['SampleID']
            allele_tbls.append(allele_tbl)

    flat = pd.concat(allele_tbls, axis='index', ignore_index=True)
    flat.to_csv(outfile, sep='\t', index=False, compression='gzip')


rule make_count_table_per_PCRrep:
    input:
        lambda wildcards:
            samplesheet.loc[samplesheet['ExperimentIDPCRRep']==wildcards.ExperimentIDPCRRep]['CRISPRessoDir']
    output:
        counts='results/byPCRRep/{ExperimentIDPCRRep}.bin_counts.txt',
        freq='results/byPCRRep/{ExperimentIDPCRRep}.bin_freq.txt'
    run:
        make_count_table(samplesheet, 'ExperimentIDPCRRep', wildcards.ExperimentIDPCRRep, get_bin_list(), output.counts, output.freq, edit_regions=config['edit_regions'])


rule make_count_table_per_experimentalRep:
    input:
        lambda wildcards:
            samplesheet.loc[samplesheet['ExperimentIDReplicates']==wildcards.ExperimentIDReplicates]['CRISPRessoDir']
    output:
        counts='results/byExperimentRep/{ExperimentIDReplicates}.bin_counts.txt',
        freq='results/byExperimentRep/{ExperimentIDReplicates}.bin_freq.txt'
    run:
        make_count_table(samplesheet, 'ExperimentIDReplicates', wildcards.ExperimentIDReplicates, get_bin_list(), output.counts, output.freq, edit_regions=config['edit_regions'])


rule trim_count_table:
    input:
        '{path}.bin_counts.txt'
    output:
        '{path}.bin_counts.topN.txt'
    params:
        n=config['max_mle_variants']
    shell:
        'head -{params.n} {input} > {output}'


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
            samplesheet.loc[samplesheet['ExperimentIDReplicates']==wildcards.ExperimentIDReplicates]['CRISPRessoDir']
    output:
        counts='results/byExperimentRepCorFilter/{ExperimentIDReplicates}.bin_counts.txt',
        freq='results/byExperimentRepCorFilter/{ExperimentIDReplicates}.bin_freq.txt'
    run:
        make_count_table(samplesheet, 'ExperimentIDReplicates', wildcards.ExperimentIDReplicates, get_bin_list(), output.counts, output.freq, input.samplesToExclude, edit_regions=config['edit_regions'])


rule make_flat_count_table_PCRrep:
    input:
        lambda wildcards: samplesheet['CRISPRessoDir']
    output:
        'results/summary/VariantCounts.flat.tsv.gz',
    run:
        make_flat_table(samplesheet, output[0])


rule make_desired_variant_tables:
    input:
        variantCounts='results/summary/VariantCounts.flat.tsv.gz',
        variantInfo=config['variant_info']
    output:
        flat="results/summary/VariantCounts.DesiredVariants.flat.tsv",
        matrix="results/summary/VariantCounts.DesiredVariants.matrix.tsv"
    params:
        codedir=config['codedir']
    shell:
        "python {params.codedir}/workflow/scripts/AggregateDesiredAlleleCounts.py --variantCounts {input.variantCounts} --variantInfo {input.variantInfo} --outbase results/summary/VariantCounts.DesiredVariants"
