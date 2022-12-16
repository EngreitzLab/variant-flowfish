#! /usr/bin/env python

import pandas as pd
import numpy as np
import argparse

def normalize_effects(raw_effects_table):
    ## Normalize effects to the most common (assumed to be wild-type) allele

    raw = raw_effects_table # pd.read_table(raw_effects_table)

    # add mean column, using formula for mean of a log-normal distribution
    raw['mean'] = np.power(10, raw['logMean']+.5*(raw['logSD']**2))
    raw['freq'] = raw['sum1'] / raw['sum1'].sum()

    # find Reference means use larger Reference mean (assuming there is both Reference and Inferred Reference) 
    refMean = raw[raw['VariantID'].str.contains('Reference')]['mean'].max()

    # normalize to reference allele
    raw['effect_size'] = raw['mean'] / refMean

    return raw

def rescale_effects_qpcr(effects_table, ff_tss_kd, qpcr_tss_kd):
    scale = qpcr_tss_kd / ff_tss_kd
    effects_table['effect_size_scaled_qpcr'] = ((effects_table['effect_size'] - 1)) * scale + 1
    effects_table.loc[effects_table['effect_size_scaled_qpcr'] < 0, 'effect_size_scaled_qpcr'] = 0 # expression can't be less than 0
    return effects_table

def adjust_effects_heterozygous_editing(effects_table, variantInfoFile, pooled='True'):
    import pdb; pdb.set_trace()
    variant_df = pd.read_table(variantInfoFile)
    num_variants = len(variant_df)
    if 'input.fraction' in effects_table.columns:
        effects_table['effect_size_scaled_adjusted'] = effects_table.apply(lambda x: 1 + get_homozygous_variant_effect_from_FF_effect(x['effect_size_scaled_qpcr'] - 1, x['input.fraction']*num_variants, pooled), axis=1)
    else:
        # if there isn't input bin, use the other bins frequency instead
        effects_table['effect_size_scaled_adjusted'] = effects_table.apply(lambda x: 1 + get_homozygous_variant_effect_from_FF_effect(x['effect_size_scaled_qpcr'] - 1, x['freq']*num_variants, pooled), axis=1)
    effects_table.loc[~(effects_table['effect_size_scaled_adjusted'] > 0), 'effect_size_scaled_adjusted'] = 0 # set effect sizes less than 0 to 0
    return effects_table

def get_homozygous_variant_effect_from_FF_effect(measured_effect, edit_rate, pooled):
    # estimate the effect size of the homozygous variant using the edit efficiency and the FF measured effect
    # this effect refers to % change, which is different from what we have in "effect size", so we subtract 1 when calling this function
    
    if pooled.lower() == 'true':
        homozygous_effect = (measured_effect * 2) / (1 + edit_rate)
    else:
        homozygous_effect = (measured_effect * 2) / (1 - (edit_rate * measured_effect))
    
    return homozygous_effect

parser = argparse.ArgumentParser()
parser.add_argument("-i", dest="input", type=str, required=True)
parser.add_argument("-o", dest="output", type=str, required=True)
parser.add_argument("-v", dest="variantInfo", type=str, required=True)
parser.add_argument("-p", dest="pooled", type=str, required=True)
parser.add_argument("-f", dest="ff_tss_guide_kd", type=float, required=True)
parser.add_argument("-q", dest="qpcr_tss_guide_kd", type=float, required=True)
args = parser.parse_args()

# read in file
raw_effects = pd.read_table(args.input)
normalized_effects = normalize_effects(raw_effects)
qpcr_rescaled_effects = rescale_effects_qpcr(normalized_effects, args.ff_tss_guide_kd, args.qpcr_tss_guide_kd)
final_effects = adjust_effects_heterozygous_editing(qpcr_rescaled_effects, args.variantInfo, args.pooled)
final_effects.to_csv(args.output, index=None, sep='\t')