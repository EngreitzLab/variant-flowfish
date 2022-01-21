import os
from random import sample
from re import L
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
        variantInfo = variantInfo[variantInfo['AmpliconID'] == sample_info['AmpliconID']] # only search variants for this amplicon
        variantSearchList = []
        
        # select rows which Aligned_Sequence contains a variant 
        unmatched_variants = []
        for index, row in variantInfo.iterrows():
            variant_df = allele_tbl[allele_tbl['Aligned_Sequence'].str.contains(row.MappingSequence)].copy()
            variant_df['Match_Sequence'] = row.MappingSequence
            variant_df['VariantID'] = row.VariantID
            variantSearchList.append(variant_df)
            # store unmatched variants
            if len(variant_df) == 0: 
                unmatched_variants.append(row)

        # combine variants and group by unique variant 
        variants = pd.concat(variantSearchList) 
        #reference_sequence = variants['Reference_Sequence'].mode().item()
        variants_grouped = variants.groupby(['Reference_Name', 'Match_Sequence', 'VariantID'])
        variant_counts = variants_grouped.sum()
        variant_counts.drop(['n_deleted', 'n_inserted', 'n_mutated'], axis=1, inplace=True)
        variant_counts['Counts'] = variants_grouped.size()
        #variant_counts['Reference_Sequence'] = reference_sequence
        variant_counts = variant_counts.reset_index()

        # get rows of allele_tbl that do not contain one of the variants
        references = allele_tbl[~allele_tbl['Aligned_Sequence'].str.contains('|'.join(variantInfo[variantInfo['RefAllele'] == False]['MappingSequence']))]
        references.drop(['n_deleted', 'n_inserted', 'n_mutated'], axis=1, inplace=True)

        # get mismatch info for references
        mismatches_column = []
        deletions_column = []
        insertions_column = []
        aligned_windows = []
        reference_windows = []
        for index, rrow in references.iterrows():
            mismatches, deletions, insertions, Aligned_Window, Reference_Window = count_mismatches(rrow['Aligned_Sequence'], rrow['Reference_Sequence'], sample_info['QuantificationWindowStart'], sample_info['QuantificationWindowEnd'])
            mismatches_column.append(mismatches)
            deletions_column.append(deletions)
            insertions_column.append(insertions)
            aligned_windows.append(Aligned_Window)
            reference_windows.append(Reference_Window)
        references['Mismatches'] = mismatches_column
        references['Deletions'] = deletions_column
        references['Insertions'] = insertions_column
        references['Aligned_Window'] = aligned_windows
        references['Reference_Window'] = reference_windows
        
        # filter out references with mismatch/insertion/deletion threshold 
        references['Errors'] = references.apply(lambda x: len(x.Mismatches)+len(x.Insertions)+len(x.Deletions), axis=1)
        reference_threshold = sample_info['ReferenceErrorThreshold']
        print('Using reference error threshold of %d' % reference_threshold)
        references = references[references['Errors'] <= reference_threshold]
                
        references.to_csv('results/variantCounts/{SampleID}.referenceAlleles.txt'.format(SampleID=sample_info['SampleID']), sep='\t', index=False)
        
        # group/sum all references together and get dict representation
        inferred_reference = references.groupby('Reference_Name')['#Reads', '%Reads'].sum().reset_index().to_dict(orient='records')[0]
        inferred_reference['VariantID'] = sample_info['AmpliconID'] + ':InferredReference'
        inferred_reference['Counts'] = len(references)
        variant_counts = variant_counts.append(inferred_reference, ignore_index=True)
        variant_counts['RefAllele'] = variant_counts['VariantID'].str.contains('Reference') # odd way to get RefAllele boolean after grouping
        
        variant_counts.rename(columns={"Reference_Name":"AmpliconID", "Match_Sequence":"MappingSequence"}, inplace=True)

        # if there are unmatched variants for this sample, add them to the variant_counts table as 0 reads
        if len(unmatched_variants) > 0:
            variant_counts = variant_counts.append(pd.DataFrame(unmatched_variants))
            variant_counts.fillna(0, inplace=True)

        variant_counts.to_csv(outfile, sep='\t', index=False)
        return

    else:  
        print('File \'%s\' not found.' % file)
        exit(1)



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
            allele_tbl = allele_tbl[['AmpliconID', 'MappingSequence', 'VariantID', 'RefAllele', '#Reads']]
            allele_tbl.rename(columns={"#Reads":row['SampleID']}, inplace=True)
            allele_tbls.append(allele_tbl)

    # combine allele_tbls into one and sum on unique variants to get count table    
    count_tbl = pd.concat(allele_tbls).groupby(['AmpliconID', 'MappingSequence', 'VariantID', 'RefAllele']).sum()
    
    # rename bins and combine duplicates 
    bin_list = bins + list(set(currSamples['Bin'].unique())-set(bins))
    combined_count_tbl = pd.DataFrame(index=count_tbl.index)
    for b in bin_list:
        rep_columns = currSamples.loc[currSamples['Bin']==b,'SampleID'].tolist()
        if (len(rep_columns) > 0):
            combined_count_tbl[b] = count_tbl[rep_columns].sum(axis=1).astype(np.int32) # sum replicate columns
        else:
            combined_count_tbl[b] = 0

    count_tbl = combined_count_tbl

    ## OUTPUT OF THIS FUNCTION:

    # one output:
    # AmpliconID    VariantID   MappingSequence RefAllele   Count   nCrispressoAlleles
    # encode the reference allele as:  HEK3 HEK3:InferredReference  [MappingSequence?]  True    [count] [nCrispressoAlleles]

    # other outputs:
    # updated count_tbl with our inferred Mismatches Deletions Insertions columns added so we can debug

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

    # extend end of current quantification window so length of window without '-' is same as reference window
    bases = RefQuantificationWindowEnd - RefQuantificationWindowStart 
    extend = 0
    i = refStart
    while bases > 0 and i < len(Reference_Sequence):
        if Reference_Sequence[i:(i+1)] != '-':
            bases -= 1
        else: 
            extend += 1
        i += 1
    refEnd = refStart + (RefQuantificationWindowEnd-RefQuantificationWindowStart) + extend

    # extend reference end if it ends in a gap
    while Reference_Sequence[refEnd-1:(refEnd)] == '-':
        refEnd = refEnd + 1
    
    start_shift = refStart - RefQuantificationWindowStart
    # end_shift = refEnd - RefQuantificationWindowEnd
        
    ## Count mismatches, deletions, and insertions in this window
    mismatches = []
    deletions = [] 
    insertions = [] 
    gap = False
    gap_spaces = 0

    for i in range(refStart,refEnd):
        if Reference_Sequence[i:(i+1)] == '-':
            if gap:
                pass
            else:
                insertions.append(i-(start_shift+gap_spaces))
                gap = True
            gap_spaces += 1 # to maintain "true" position of quantification window

        elif Aligned_Sequence[i:(i+1)] == '-':
            if gap:
                pass
            else:
                deletions.append(i-(start_shift+gap_spaces))
                gap = True
                # gap_spaces += 1

        elif Reference_Sequence[i:(i+1)] != Aligned_Sequence[i:(i+1)]:
            mismatches.append(i-(start_shift+gap_spaces))
            gap = False

        else:
            gap = False

    Aligned_Window = Aligned_Sequence[refStart:refEnd]
    Reference_Window = Reference_Sequence[refStart:refEnd]
    
    return (mismatches, deletions, insertions, Aligned_Window, Reference_Window)

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
    parser.add_argument("--reference_threshold", required=False, default=None)
    parser.add_argument("--edit_regions", required=False, default=None)

    args = parser.parse_args()
    aggregate_variant_counts(args.samplesheet, args.sample_id, args.outfile, args.variantInfoFile, args.reference_threshold)
    make_count_table(args.samplesheet, args.group_col, args.group_id, args.bins, args.outfile, args.outfile_frequencies, args.samplesToExclude, args.edit_regions)
    
if __name__ == "__main__":
    main()
