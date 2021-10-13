import argparse
import os
import sys
import pandas as pd


def getAlleleTable(countsFlat, variantInfo, minFreq):
    '''
    Takes in DataFrames of variant count flat file and desired variant data and
    performs matching of variant count data with desired variants, combining
    read data for sequences matching single variants.
    Returns flat and matrix versions of desired counts tables with combined read
    info, as well as data for the aligned sequences that matched multiple variants.
    '''
    
    # iterate through variants from variantInfo and add matches from
    # countsFlat to dataframe list
    variantSearchList = []
    count = 0
    for variant in variantInfo.MappingSequence:
        count += 1
        if count % 100 == 0:
            print(count, "variants processed")
        variant_df = countsFlat[countsFlat['Aligned_Sequence'].str.contains(
            variant)]
        variant_df['MatchSequence'] = variant
        variantSearchList.append(variant_df)

    variantMatches = pd.concat(variantSearchList)

    # group matches by consistent columns and list #Reads, %Reads,
    # MatchSequence to see if there are multiple matches
    variantsGrouped = variantMatches.groupby(
        [
            'Aligned_Sequence',
            'Reference_Sequence',
            'n_deleted',
            'n_inserted',
            'n_mutated',
            'SampleID'],
        as_index=False,
        observed=True)[
            '#Reads',
            '%Reads',
            'MatchSequence'].agg(
                lambda x: list(x))

    print('number of variant matches:')
    print(variantsGrouped['MatchSequence'].map(len).value_counts())
    matchCounts = pd.DataFrame(
        variantsGrouped['MatchSequence'].map(len).value_counts()).reset_index()
    matchCounts.columns = ['num_matches', 'count']

    # collect multiple matches
    variantsGroupedMultiple = variantsGrouped[variantsGrouped['MatchSequence'].map(
        len) > 1]
    variantsGroupedMultiple['num_matches'] = variantsGroupedMultiple['MatchSequence'].apply(
        lambda x: len(x))

    # keep single variant matches
    variantsGroupedSingle = variantsGrouped[variantsGrouped['MatchSequence'].map(
        len) == 1]
    variantsGroupedSingle['#Reads'] = variantsGroupedSingle['#Reads'].apply(
        lambda x: x[0])
    variantsGroupedSingle['%Reads'] = variantsGroupedSingle['%Reads'].apply(
        lambda x: x[0])
    variantsGroupedSingle['MatchSequence'] = variantsGroupedSingle['MatchSequence'].apply(
        lambda x: x[0])

    # group by Sample ID and MatchSequence, sum up the # and % reads, reset
    # DataFrame index, merge with all variants
    desired_counts_flat = variantsGroupedSingle.groupby(
        ['SampleID', 'MatchSequence'])[['#Reads', '%Reads']].sum().reset_index()
    desired_counts_flat = desired_counts_flat.merge(
        variantInfo,
        left_on='MatchSequence',
        right_on='MappingSequence',
        how='outer').fillna(0)
    desired_counts_flat = desired_counts_flat[['MappingSequence',
                                               'MatchSequence',
                                               'AmpliconID',
                                               'VariantID',
                                               'RefAllele',
                                               'SampleID',
                                               '%Reads',
                                               '#Reads']]

    desired_counts_table = desired_counts_flat.pivot(
        index=[
            'MatchSequence',
            'AmpliconID',
            'VariantID',
            'RefAllele'],
        columns='SampleID',
        values='#Reads').reset_index()
    desired_counts_table = desired_counts_table.rename_axis(None, axis=1).reset_index(
        drop=True).fillna(0)  # the index needs to be reset because pivot names the index "SampleID"

    # filter on minFrequency
    desired_counts_table = desired_counts_table[desired_counts_table.max(
        axis=1) > minFreq]

    return(desired_counts_flat, desired_counts_table, variantsGroupedMultiple, matchCounts)


def main():
    parser = argparse.ArgumentParser(
        description='Aggregate allele frequency information across CRISPResso runs into a table for plotting and analysis.')
    parser.add_argument(
        '--variantCounts',
        type=str,
        required=True,
        help='Variant count flat file (gzipped) from the variant-flowfish snakemake pipeline')
    parser.add_argument(
        '--variantInfo',
        type=str,
        required=True,
        help='File containing desired variants. Tab-delimited file containing columns AmpliconID, GuideSpacer, MappingSequence, RefAllele;  where MappingSequence matches the Aligned_Sequence output column in the Alleles frequency table in CRISPResso')
    parser.add_argument(
        '--outbase',
        type=str,
        required=True,
        help='Output filebase of allele frequencies (alleles x samples)')
    parser.add_argument(
        '--minFrequency',
        type=int,
        default=0,
        help='Minimum threshold on alleles to include in merged table')

    args = parser.parse_args()

    if not os.path.exists(args.variantCounts):
        sys.exit('AggregateAlleleCounts: --variantCounts file not found.')

    if not os.path.exists(args.variantInfo):
        sys.exit('AggregateAlleleCounts: --variantInfo file not found.')

    countsFlat = pd.read_csv(args.variantCounts, sep='\t')
    variantInfo = pd.read_csv(args.variantInfo, sep='\t')

    desiredCounts = getAlleleTable(countsFlat, variantInfo, minFreq=args.minFrequency)

    desiredCounts[0].to_csv(args.outbase + '.flat.tsv', sep='\t', index=False)
    desiredCounts[1].to_csv(args.outbase + '.matrix.tsv', sep='\t', index=False)
    desiredCounts[2].to_csv(args.outbase + '.multiples.tsv', sep='\t', index=False)
    desiredCounts[3].to_csv(args.outbase + '.match_counts.tsv', sep='\t', index=False)


if __name__ == "__main__":
    main()
