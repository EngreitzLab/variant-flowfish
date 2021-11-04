#! /usr/bin/env python

import pandas as pd
import numpy as np
import argparse

def normalize_effects(raw_effects_file, outfile, variantInfo=None, index="MappingSequence"):
    ## Normalize effects to the most common (assumed to be wild-type) allele

    # read in files
    raw = pd.read_table(raw_effects_file)

    # add mean column, using formula for mean of a log-normal distribution
    raw['mean'] = np.power(10, raw['logMean']+.5*(raw['logSD']**2))
    raw['freq'] = raw['sum1'] / raw['sum1'].sum()

    if variantInfo is None:
        # get reference allele expressions, assuming the most frequent allele is the reference allele
        refMean = raw.sort_values('sum1', axis=1, ascending=False)['mean'].values[0]
    else:
        # merge in variant info to raw effects table 
        variants = pd.read_table(variantInfo).set_index(index)
        RefAllele_variants = variants[variants['RefAllele'] == True].reset_index()

       # search for variants as substrings of raw effects MappingSequence
        variantSearchList = []
        count = 0
        for index, row in RefAllele_variants.iterrows():
            count += 1
            if count % 100 == 0:
                print(count, "variants processed")
            variant_df = raw[raw['MappingSequence'].str.contains(row.MappingSequence)]
            variant_df['MatchSequence'] = row.MappingSequence
            variant_df['VariantID'] = row.VariantID
            variantSearchList.append(variant_df)
        variantMatches = pd.concat(variantSearchList)
        
        # normalize to most abundant annotated reference allele
        refMean = variantMatches.sort_values('sum1', axis=0, ascending=False)['mean'].values[0]
        raw = variantMatches # not sure but seems like need to keep only variant matches

    # normalize to reference allele
    raw['effect_size'] = raw['mean'] / refMean

    # write out
    raw.to_csv(outfile, index=None, sep='\t')


parser = argparse.ArgumentParser()
parser.add_argument("-i", dest="input", type=str, required=True)
parser.add_argument("-o", dest="output", type=str, required=True)
parser.add_argument("-v", dest="variantInfo", type=str, required=False, default=None)
args = parser.parse_args()

normalize_effects(args.input, args.output, args.variantInfo)
