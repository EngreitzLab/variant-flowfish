#! /usr/bin/env python

import pandas as pd
import numpy as np
import argparse

def normalize_effects(raw_effects_file, outfile, variantInfo=None, index="MappingSequence"):
    ## Normalize effects to the most common (assumed to be wild-type) allele

    # read in files
    raw = pd.read_table(raw_effects_file).set_index(index)

    if variantInfo is not None:
        variants = pd.read_table(variantInfo).set_index(index)
        raw = raw.join(variants, how='left', on=index)

    # add mean column, using formula for mean of a log-normal distribution
    raw['mean'] = np.power(10, raw['logMean']+.5*(raw['logSD']**2))
    raw['freq'] = raw['sum1'] / raw['sum1'].sum()
    if variantInfo is None:
        # get reference allele expressions, assuming the most frequent allele is the reference allele
        refMean = raw.sort_values('sum1', axis=1, ascending=False)['mean'].values[0]
    else:
        # normalize to most abundant annotated reference allele
        refMean = raw.loc[raw['RefAllele'].values == True,].sort_values('sum1', axis=0, ascending=False)['mean'].values[0]

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
