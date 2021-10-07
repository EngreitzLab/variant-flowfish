import argparse
import os
import sys
import pandas as pd
import warnings


def getAlleleTable(countsFlat, variantInfo, minFreq):
    countsFiltered = df.copy(countsFlat)

    # find all instances of VariantInfo mapping sequences in each Aligned_Sequence
    countsFiltered['Matches'] = countsFiltered['Aligned_Sequence'].str.findall('|'.join(variantInfo.MappingSequence))

    # check if there is one match for each or not
    if not all(countsFiltered['Matches'].str.len() == 1):
        warnings.warn('AggregateAlleleCounts: Sequence exists with either none or multiple matches.')
    # countsFlat[countsFlat['MatchSequence'].str.len() != 1] # if you want to get the entries with none/multiple matches

    # keep first result from Matches
    countsFlat['MatchSequence'] = countsFlat['Matches'].str[0]

    # group by Sample ID and MatchSequence, sum up the # and % reads, reset DataFrame index
    desired_counts_flat = countsFiltered.groupby(['SampleID', 'MatchSequence'])[['#Reads', '%Reads']].sum().reset_index()

    # create table with SampleIDs and #Reads
    desired_counts_table = desired_counts_flat.pivot(index='MatchSequence',columns='SampleID',values='#Reads').reset_index()
    desired_counts_table = desired_counts_table.rename_axis(None, axis=1).reset_index(drop=True) # the index needs to be reset because pivot names the index "SampleID"


def main():
    parser = argparse.ArgumentParser(description='Aggregate allele frequency information across CRISPResso runs into a table for plotting and analysis.')
    parser.add_argument('--variantCounts', type=str, required=True,
                        help='Variant count flat file (gzipped) from the variant-flowfish snakemake pipeline')
    parser.add_argument('--variantInfo', type=str, required=True,
                        help='File containing desired variants. Tab-delimited file containing columns AmpliconID, GuideSpacer, MappingSequence, RefAllele;  where MappingSequence matches the Aligned_Sequence output column in the Alleles frequency table in CRISPResso')
    parser.add_argument('--outbase', type=str, required=True,
                        help='Output filebase of allele frequencies (alleles x samples)')
    parser.add_argument('--minFrequency', type=int, default=0,
                        help='Minimum threshold on alleles to include in merged table')

    args = parser.parse_args()

    if not os.path.exists(args.variantCounts):
        sys.exit('AggregateAlleleCounts: --variantCounts file not found.')

    if not os.path.exists(args.variantInfo):
        sys.exit('AggregateAlleleCounts: --variantInfo file not found.')

    countsFlat = pd.read_csv(args.variantCounts, sep='\t')
    variantInfo = pd.read_csv(args.variantInfo, sep='\t')

    desiredCounts = getAlleleTable(countsFlat, variantInfo, minFreq=args.minFrequency)

if __name__=="__main__":
    main()
