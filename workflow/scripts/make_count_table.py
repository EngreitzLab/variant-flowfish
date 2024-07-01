import os
from random import sample
from re import L
import pandas as pd
import numpy as np
import argparse
import regex as re
import ahocorasick

pd.options.mode.chained_assignment = None  # default='warn', turn off warnings bc we are doing assignments correctly

# make table of all variants
def make_flat_table(samplesheet, outfile):
    allele_tbls = []
    for idx, row in samplesheet.iterrows():
        file = row['variantCountFile']

        if (os.path.exists(file)):
            try:
                allele_tbl = pd.read_csv(file, sep='\t')
                allele_tbl['SampleID'] = row['SampleID']
                allele_tbl['BioRep'] = row['BioRep']

                if 'FlowFISHRep' in row.keys():
                    allele_tbl['FlowFISHRep'] = row['FlowFISHRep']
                else:
                    print('No FlowFISHRep column (ok if genotyping)')

                allele_tbl['PCRRep'] = row['PCRRep']

                if 'Bin' in row.keys():
                    allele_tbl['Bin'] = row['Bin']
                else:
                    print('No Bin column (ok if genotyping)')

                # splitting like this is specific to our current syntax for VariantID, may need to change
                allele_tbl['Location'] = allele_tbl['VariantID'].str.split(':').str[:3].str.join(':')
                allele_tbl['Variant'] = allele_tbl['VariantID'].str.split(':').str[-1:].str.join(':')
                allele_tbls.append(allele_tbl)
            except:
                print("Error reading file: " + file)
                continue
        else:
            print('File not found:', file)

    try:
        flat = pd.concat(allele_tbls, axis='index', ignore_index=True)
    except:
        # need to error handle this?
        flat = pd.DataFrame()
    flat.to_csv(outfile, sep='\t', index=False, compression='gzip')

# returns list of indexes where the variant sequence was found
def find_variant_locations(variant_sequence, aligned_sequence, reference_sequence):
    # potentially update to accommodate for gaps and indels
    locs = np.array([i.start() for i in re.finditer(variant_sequence, aligned_sequence)])
    for i in range(len(locs)):
        num_gaps = reference_sequence[:locs[i]].count('-')
        locs[i] -= num_gaps
    return list(set(locs))

def aggregate_variant_counts(samplesheet, SampleID, outfile, variantInfoFile):
    sample_info = samplesheet.loc[samplesheet['SampleID'] == SampleID].to_dict('records')[0]
    file = "results/crispresso/CRISPResso_on_{SampleID}/Alleles_frequency_table.zip".format(SampleID=SampleID)
    if (os.path.exists(file)):
        allele_tbl = pd.read_csv(file, sep='\t')
        allele_tbl['#Reads'] = allele_tbl['#Reads'].astype(np.int32)
        variantInfo = pd.read_table(variantInfoFile)
        variantInfo = variantInfo[['AmpliconID', 'VariantID', 'MappingSequence', 'RefAllele']] # keep required cols only
        variantInfo = variantInfo[variantInfo['AmpliconID'] == sample_info['AmpliconID']] # only search variants for this amplicon
        if len(variantInfo) == 0:
            print('No matching variants found for AmpliconID: %s, check inputs.' % sample_info['AmpliconID'])
            exit(0)
        variantSearchList = []        

        # Create an Aho-Corasick Automaton
        A = ahocorasick.Automaton()

        # Make MappingSequence - VariantID dictionary
        variant_dict = pd.DataFrame(variantInfo[['MappingSequence', 'VariantID']]).set_index('MappingSequence')['VariantID'].to_dict() 
        for idx, variant in enumerate(variant_dict.keys()):
            A.add_word(variant, (idx, variant))
        A.make_automaton()

        matches = []
        match_counts = []
        # Iterate over all reads
        for index, row in allele_tbl.iterrows():
            num_matches = 0
            # Search for
            results = A.iter(row['Aligned_Sequence'])        

            # Iterate over the results
            variant_list = []
            for end_index, (idx, word) in results:
                num_matches += 1
                variant_list.append(variant_dict[word]) 
            if num_matches > 1:
                matches.append(variant_list)
            elif num_matches == 1:
                matches.append(variant_list[0])
            else: # no match 
                matches.append(None)
            match_counts.append(num_matches)
        allele_tbl['n_variant_matches'] = match_counts
        allele_tbl['VariantID'] = matches

        # allele_tbl_count_data = pd.DataFrame(allele_tbl.groupby(['num_variant_matches'])['#Reads'].sum()).T
        # SampleID_col = [SampleID]*len(allele_tbl)
        # allele_tbl.insert(loc=0, column='SampleID', value=SampleID_col) 
        # allele_tbl.to_csv('read_data.csv', mode='a', index=False, header=False)

        unmatched_variants = []

        for index, row in variantInfo.iterrows():
            variant_df = allele_tbl[allele_tbl['VariantID'] == row['VariantID']]
            # If allele_tbl['VariantID'] exactly matches a VariantID, it is the only match in that read

            # store unmatched variants and continue to next variant in iteration
            if len(variant_df) == 0:
                unmatched_variants.append(row)
                continue

            variant_df['Match_Sequence'] = row.MappingSequence
            variantSearchList.append(variant_df)

        if len(variantSearchList) > 0:
            # combine variants and group by unique variant
            variants = pd.concat(variantSearchList)
            variants_grouped = variants.groupby(['Reference_Name', 'Match_Sequence', 'VariantID'])
            variant_counts = variants_grouped.sum()
            variant_counts.drop(['n_deleted', 'n_inserted', 'n_mutated'], axis=1, inplace=True)
            variant_counts['Counts'] = variants_grouped.size()
            variant_counts = variant_counts.reset_index()
        else:
            # set up empty variant_counts dataframe if no variants found
            variant_counts = pd.DataFrame(columns=['Reference_Name', 'Match_Sequence', 'VariantID', '#Reads', '%Reads', 'Counts', 'MatchLocations', 'RefAllele'])

        # get rows of allele_tbl that do not contain one of the variants
        references = allele_tbl[allele_tbl['n_variant_matches'] == 0]
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

        if len(references) == 0: # no reference alleles found
            references['Errors'] = []

        else:
            # filter out references with mismatch/insertion/deletion threshold
            references['Errors'] = references.apply(lambda x: len(x.Mismatches)+len(x.Insertions)+len(x.Deletions), axis=1)

        reference_threshold = sample_info['ReferenceErrorThreshold']
        print('Using reference error threshold of %d' % reference_threshold)

        # save all references before filtering
        references.to_csv('results/variantCounts/{SampleID}.referenceAlleles.txt'.format(SampleID=sample_info['SampleID']), sep='\t', index=False)

        references = references[references['Errors'] <= reference_threshold]

        # group/sum all references together and get dict representation
        grouped_references = references.groupby('Reference_Name')[['#Reads', '%Reads']].sum().reset_index().to_dict(orient='records')
        if len(grouped_references) == 1: # expected behavior, all references should compress to one entry
            inferred_reference = grouped_references[0]
            inferred_reference['VariantID'] = sample_info['AmpliconID'] + ':InferredReference'
            inferred_reference['Counts'] = len(references)
        else:
            print('No references found after filtering.')
            inferred_reference = pd.DataFrame(columns=['Reference_Name', 'Match_Sequence', 'VariantID', '#Reads', '%Reads', 'Counts', 'MatchLocations', 'RefAllele'])
            inferred_reference['VariantID'] = sample_info['AmpliconID'] + ':InferredReference'
            inferred_reference['Counts'] = 0
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
        print('File \'%s\' not found; writing empty variant and reference files.' % file)
        empty_variant = pd.DataFrame(columns=['AmpliconID', 'MappingSequence', 'VariantID', '#Reads', '%Reads', 'Counts', 'MatchLocations', 'RefAllele'])
        empty_variant.to_csv(outfile, sep='\t', index=False)
        empty_reference = pd.DataFrame(columns=['Aligned_Sequence', 'Reference_Sequence', 'Reference_Name', 'Read_Status', '#Reads', '%Reads', 'Mismatches', 'Deletions', 'Insertions', 'Aligned_Window', 'Reference_Window', 'Errors'])
        empty_reference.to_csv('results/variantCounts/{SampleID}.referenceAlleles.txt'.format(SampleID=sample_info['SampleID']), sep='\t', index=False)
        return
        # exit(1)



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
            try:
                allele_tbl = pd.read_csv(file, sep='\t')
            except:
                # empty file
                print('Error reading file: %s' % file)
                continue

            if len(allele_tbl) == 0: # empty table
                print('%s is empty.' % file)
                continue

            allele_tbl['#Reads'] = allele_tbl['#Reads'].astype(np.int32)

            # pare down allele_tbl columns for counts
            allele_tbl = allele_tbl[['AmpliconID', 'MappingSequence', 'VariantID', 'RefAllele', '#Reads']]
            allele_tbl.rename(columns={"#Reads":row['SampleID']}, inplace=True)
            allele_tbls.append(allele_tbl)

    try:
        # combine allele_tbls into one and sum on unique variants to get count table
        # fill NaNs to make sure nothing dropped from grouping
        all_allele_tbls = pd.concat(allele_tbls).fillna(0)
        all_allele_tbls['MappingSequence'] = all_allele_tbls['MappingSequence'].astype(str) # handle string zero vs int zero
        count_tbl = all_allele_tbls.groupby(['AmpliconID', 'MappingSequence', 'VariantID', 'RefAllele']).sum()

        # rename bins and combine duplicates
        bin_list = bins + list(set(currSamples['Bin'].unique())-set(bins))
        combined_count_tbl = pd.DataFrame(index=count_tbl.index)
        for b in bin_list:
            rep_columns = currSamples.loc[currSamples['Bin']==b,'SampleID'].tolist()
            if (len(rep_columns) > 0):
                combined_count_tbl[b] = count_tbl[count_tbl.columns.intersection(rep_columns)].sum(axis=1).astype(np.int32) # sum replicate columns, only selecting ones in the count tbl (could be that variant counts file was empty, so column isn't there)
            else:
                combined_count_tbl[b] = 0

        count_tbl = combined_count_tbl
    except:
        print("Error combining allele tables, writing empty count table.")
        count_tbl = pd.DataFrame(columns=['AmpliconID', 'MappingSequence', 'VariantID', 'RefAllele'])

    ## OUTPUT OF THIS FUNCTION:

    # one output:
    # AmpliconID    VariantID   MappingSequence RefAllele   Count   nCrispressoAlleles
    # encode the reference allele as:  HEK3 HEK3:InferredReference  [MappingSequence?]  True    [count] [nCrispressoAlleles]

    # other outputs:
    # updated count_tbl with our inferred Mismatches Deletions Insertions columns added so we can debug

    count_tbl.to_csv(outfile, sep='\t')

    freq_tbl = count_tbl.reset_index()

    if len(freq_tbl[freq_tbl['RefAllele'] == True]) > 1: # avoid double counting "true" and "inferred" reference in freq table
        freq_tmp_no_infer = freq_tbl[~freq_tbl['VariantID'].str.contains('InferredReference')][count_tbl.columns] # select count data without inferred reference
        freq_tmp = freq_tbl[count_tbl.columns] # save count data with inferred 
        freq_tbl[count_tbl.columns] = freq_tmp.div(freq_tmp_no_infer.sum(axis=0), axis=1) # freq divide; not including Inferred Reference count in sum
        freq_tbl.set_index(count_tbl.index.names, inplace=True)  
    else:
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

    ## Count mismatches, deletions, and insertions in this window
    mismatches = []
    deletions = []
    insertions = []
    gap = False
    gap_spaces = 0

    for i in range(refStart,refEnd):
        # added lines to count every base insertion and deletion, not just the first instance
        if Reference_Sequence[i:(i+1)] == '-':
            if gap:
                insertions.append(i-(start_shift+gap_spaces))
                pass
            else:
                insertions.append(i-(start_shift+gap_spaces))
                gap = True
            gap_spaces += 1 # to maintain "true" position of quantification window

        elif Aligned_Sequence[i:(i+1)] == '-':
            if gap:
                deletions.append(i-(start_shift+gap_spaces))
                pass
            else:
                deletions.append(i-(start_shift+gap_spaces))
                gap = True

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
    parser.add_argument("--edit_regions", required=False, default=None)

    args = parser.parse_args()
    aggregate_variant_counts(args.samplesheet, args.sample_id, args.outfile, args.variantInfoFile)
    make_count_table(args.samplesheet, args.group_col, args.group_id, args.bins, args.outfile, args.outfile_frequencies, args.samplesToExclude, args.edit_regions)

if __name__ == "__main__":
    main()
