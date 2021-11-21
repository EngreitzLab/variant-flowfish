import os
import pandas as pd
import numpy as np
import argparse

def make_count_table(samplesheet, group_col, group_id, bins, outfile, outfile_frequencies, variantInfo=None, samplesToExclude=None, edit_regions=None):
    ## Function to make a count table at various layers of resolution (e.g., by experiment, or by replicate, or by PCR replicate)

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
        # count_tbl = count_tbl.set_index('Aligned_Sequence')
    
    # merge with variants table (maybe make separate function?)
    count_tbl = count_tbl.reset_index() # resetting index to use Aligned_Sequence column
    variants = pd.read_table(variantInfo)
    variantSearchList = []
    for index, row in variants.iterrows():
        variant_df = count_tbl[count_tbl['Aligned_Sequence'].str.contains(row.MappingSequence)]
        variant_df['MatchSequence'] = row.MappingSequence
        variant_df['VariantID'] = row.VariantID
        variantSearchList.append(variant_df)

    # make count_tbl the variant matches list and group by unique variant 
    count_tbl = pd.concat(variantSearchList)
    count_tbl = count_tbl.groupby(['MatchSequence', 'VariantID']).sum()
    count_tbl = count_tbl.reset_index()
    count_tbl.rename(columns={'MatchSequence':'MappingSequence'}, inplace=True)

    count_tbl.to_csv(outfile, sep='\t', index=False)
    
    count_tbl = count_tbl.set_index(['MappingSequence', 'VariantID']) # needed for freq operation
    freq_tbl = count_tbl.div(count_tbl.sum(axis=0), axis=1)
    freq_tbl.to_csv(outfile_frequencies, sep='\t', float_format='%.6f')



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--group_col", type=str, required=True)
    parser.add_argument("--group_id", type=str, required=True) 
    parser.add_argument("--bins", required=True) 
    parser.add_argument("--outfile", type=str, required=True) 
    parser.add_argument("--outfile_frequencies", type=str, required=True) 
    parser.add_argument("--samplesToExclude", required=False, default=None)
    parser.add_argument("--edit_regions", required=False, default=None)

    args = parser.parse_args()
    make_count_table(args.samplesheet, args.group_col, args.group_id, args.bins, args.outfile, args.outfile_frequencies, args.samplesToExclude, args.edit_regions)

if __name__ == "__main__":
    main()
