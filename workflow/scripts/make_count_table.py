import os
from random import sample
import pandas as pd
import numpy as np
import argparse
import regex as re

def aggregate_variant_counts(samplesheet, SampleID, outfile, variantInfoFile):
    sample_info = samplesheet.loc[samplesheet['SampleID'] == SampleID].to_dict('records')[0] 
    file = "results/crispresso/CRISPResso_on_{SampleID}/Alleles_frequency_table.zip".format(SampleID=SampleID)
    if (os.path.exists(file)):
        allele_tbl = pd.read_csv(file, sep='\t')
        allele_tbl['#Reads'] = allele_tbl['#Reads'].astype(np.int32)
        variantInfo = pd.read_table(variantInfoFile)
        variantSearchList = []
        
        # select rows which Aligned_Sequence contains a variant 
        for index, row in variantInfo.iterrows():
            variant_df = allele_tbl[allele_tbl['Aligned_Sequence'].str.contains(row.MappingSequence)]
            variant_df['Match_Sequence'] = row.MappingSequence
            variant_df['VariantID'] = row.VariantID
            variantSearchList.append(variant_df)
        
        # combine variants and group by unique variant 
        variants = pd.concat(variantSearchList) 
        reference_sequence = variants['Reference_Sequence'].mode().item()
        variants_grouped = variants.groupby(['Reference_Name', 'Match_Sequence', 'VariantID'])
        variant_counts = variants_grouped.sum()
        variant_counts.drop(['n_deleted', 'n_inserted', 'n_mutated'], axis=1, inplace=True)
        variant_counts['Counts'] = variants_grouped.size()
        variant_counts['Reference_Sequence'] = reference_sequence
        variant_counts = variant_counts.reset_index()

        # get rows of allele_tbl that do not contain one of the variants
        references = allele_tbl[~allele_tbl['Aligned_Sequence'].str.contains('|'.join(variantInfo[variantInfo['RefAllele'] == False]['MappingSequence']))]
        
        # get mismatch info for references
        nMismatches_column = []
        nDeletions_column = []
        nInsertions_column = []
        for index, rrow in references.iterrows():
            nMismatches, nDeletions, nInsertions = count_mismatches(rrow['Aligned_Sequence'], rrow['Reference_Sequence'], sample_info['QuantificationWindowStart'], sample_info['QuantificationWindowEnd'])
            nMismatches_column.append(nMismatches)
            nDeletions_column.append(nDeletions)
            nInsertions_column.append(nInsertions)
        references['nMismatches'] = nMismatches_column
        references['nDeletions'] = nDeletions_column
        references['nInsertions'] = nInsertions_column
        references.to_csv('results/variantCounts/{SampleID}.referenceAlleles.txt'.format(SampleID=sample_info['SampleID']), sep='\t', index=False)

        # group/sum all references together and get dict representation
        inferred_reference = references.groupby('Reference_Name')['#Reads', '%Reads'].sum().reset_index().to_dict(orient='records')[0]
        inferred_reference['VariantID'] = sample_info['AmpliconID'] + ':InferredReference'
        inferred_reference['Counts'] = len(references)
        inferred_reference['Reference_Sequence'] = references['Reference_Sequence'].mode().item()
        variant_counts = variant_counts.append(inferred_reference, ignore_index=True)
        variant_counts['RefAllele'] = variant_counts['VariantID'].str.contains('Reference') # odd way to get RefAllele boolean after grouping

    variant_counts.to_csv(outfile, sep='\t', index=False)

def make_count_table(samplesheet, group_col, group_id, bins, outfile, outfile_frequencies, variantInfo=None, samplesToExclude=None):
    ## Function to make a count table at various layers of resolution (e.g., by experiment, or by replicate, or by PCR replicate)
    currSamples = samplesheet.loc[samplesheet[group_col]==group_id]
    currSamples['Bin'] = currSamples['Bin'].fillna('NA')
    if samplesToExclude is not None:
        exclude = pd.read_table(samplesToExclude, header=None, names=['SampleID'])
        currSamples = currSamples[~currSamples['SampleID'].isin(exclude['SampleID'].values)]

    allele_tbls = []
    for idx, row in currSamples.iterrows():
        file = "results/variantCounts/{SampleID}.variantCounts.txt".format(SampleID=row['SampleID'])
        if os.path.exists(file):
            allele_tbl = pd.read_csv(file, sep='\t')
            allele_tbl['#Reads'] = allele_tbl['#Reads'].astype(np.int32)

            # pare down allele_tbl columns for counts
            allele_tbl = allele_tbl[['Reference_Name', 'Match_Sequence', 'VariantID', 'RefAllele', '#Reads']]
            allele_tbl.rename(columns={"#Reads":row['SampleID']}, inplace=True)
            allele_tbls.append(allele_tbl)

    # combine allele_tbls into one and sum on unique variants to get count table    
    count_tbl = pd.concat(allele_tbls).groupby(['Reference_Name', 'Match_Sequence', 'VariantID', 'RefAllele']).sum()
    
    # rename bins and combine duplicates 
    bin_list = bins + list(set(currSamples['Bin'].unique())-set(bins))
    combined_count_tbl = pd.DataFrame(index=count_tbl.index)
    for b in bin_list:
        regex_bin = '.*Bin{}$'.format(b) # make regex to select column ending in BinA, BinB, etc.
        rep_columns = count_tbl.filter(regex=regex_bin).columns
        combined_count_tbl[b] = count_tbl[rep_columns].sum(axis=1) # sum replicate columns
    count_tbl = combined_count_tbl

    ## OUTPUT OF THIS FUNCTION:

    # one output:
    # AmpliconID    VariantID   MappingSequence RefAllele   Count   nCrispressoAlleles
    # encode the reference allele as:  HEK3 HEK3:InferredReference  [MappingSequence?]  True    [count] [nCrispressoAlleles]

    # other outputs:
    # updated count_tbl with our inferred nMismatches nDeletions nInsertions columns added so we can debug

    count_tbl.to_csv(outfile, sep='\t')
    
    freq_tbl = count_tbl.div(count_tbl.sum(axis=0), axis=1)
    freq_tbl.to_csv(outfile_frequencies, sep='\t', float_format='%.6f')

def count_mismatches(Aligned_Sequence, Reference_Sequence, RefQuantificationWindowStart, RefQuantificationWindowEnd):
    ## this function is a bit tricky, has to account for gaps/indels in the alignment marked by '-'
    ## I wrote some code here but this needs to be validated
    ## Find the coordinates to analyze, accounting for gaps/indels in the alignment marked by '-'
    refGaps = [m.start() for m in re.finditer('-', Reference_Sequence)]
    numGapsBeforeQws = (pd.Series(refGaps) < RefQuantificationWindowStart).sum()
    refStart = RefQuantificationWindowStart + numGapsBeforeQws
    
    while Reference_Sequence[refStart:(refStart+1)] == '-':
        refStart = refStart + 1

    refEnd = refStart + (RefQuantificationWindowEnd-RefQuantificationWindowStart)

    while Reference_Sequence[refEnd:(refEnd+1)] == '-':
        refEnd = refEnd + 1
    
    ## Count mismatches, deletions, and insertions in this window
    nMismatches=0
    nDeletions=0
    nInsertions=0
    for i in range(refStart,refEnd):
        if Reference_Sequence[i:(i+1)] == '-':
            nInsertions = nInsertions + 1
        elif Aligned_Sequence[i:(i+1)] == '-':
            nDeletions = nDeletions + 1
        elif Reference_Sequence[i:(i+1)] != Aligned_Sequence[i:(i+1)]:
            nMismatches = nMismatches + 1
    
    return (nMismatches, nDeletions, nInsertions)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--group_col", type=str, required=True)
    parser.add_argument("--group_id", type=str, required=True) 
    parser.add_argument("--bins", required=True) 
    parser.add_argument("--outfile", type=str, required=True) 
    parser.add_argument("--outfile_frequencies", type=str, required=True) 
    parser.add_argument("--variantInfoFile", type=str, required=True)
    parser.add_argument("--samplesToExclude", required=False, default=None)
    parser.add_argument("--edit_regions", required=False, default=None)

    args = parser.parse_args()
    aggregate_variant_counts(args.samplesheet, args.sample_id, args.outfile, args.variantInfoFile)
    make_count_table(args.samplesheet, args.group_col, args.group_id, args.bins, args.outfile, args.outfile_frequencies, args.samplesToExclude, args.edit_regions)

if __name__ == "__main__":
    main()
