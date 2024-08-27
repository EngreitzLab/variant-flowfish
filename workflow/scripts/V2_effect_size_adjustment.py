#! /usr/bin/env python

import pandas as pd
import numpy as np
import argparse

def adjust_effects_heterozygous_editing_V2(effects_table, output_file):
	
	#read in V1 effect size table
	effects_table = pd.read_table(effects_table, compression='gzip')	

	# remove WT alleles, we dont need them for this 
	effects_table = effects_table[effects_table['RefAllele'] != True]
	# set threshold for frequency of edits in the population that we trust effect sizes for
	effects_table = effects_table[effects_table['freq'] >= 0.0001] 

	# set variables
	loops = 0 # counter
	thres = 0.001 # threshold for expectation maximization. Once met it will break the loop
	effects_table.rename(columns={'effect_size': 'V1_effect_size'}, inplace=True) # rename the initial MLE effect size column so it doesnt get overwritten
	effects_table['effect_size_updated'] = effects_table['effect_size_scaled_qpcr'] -1 # ni in the space of the z0, the measured effect from original MLE scaled for qPCR, or e_new
	effects_table['effect_size_eold'] = -0.20 # random number, here just starting with 0. this is the e_old
	while ((effects_table['effect_size_eold'] - effects_table['effect_size_updated'])**2).mean() > thres:
		effects_table['effect_size_eold'] = effects_table['effect_size_updated'] # setting for each loop iteration, e_old=e_new

		# Step 1: Calculate WT_drag_contribution and ref_allele_contribution for each row
		effects_table['WT_drag_contribution'] = (effects_table['guide_freq'] * (1 - effects_table['adjusted_freq'])) * (1 + (effects_table['effect_size_eold'] * effects_table['adjusted_freq']) / 2) #will sum these in next step
		effects_table['ref_allele_contribution'] = effects_table['guide_freq'] * (1 - effects_table['adjusted_freq']) #will sum these in next step

		# Step 2: Calculate WT_drag_contribution_summed, f_0, and d for each sample without using groupby
		unique_combinations = effects_table[['BioRep', 'FFRep']].drop_duplicates()
		for index, row in unique_combinations.iterrows():
			BioRep = row['BioRep']
			FFRep = row['FFRep']
			# Create subset DataFrame for the current BioRep and FFRep combination
			subset_df = effects_table[(effects_table['BioRep'] == BioRep) & (effects_table['FFRep'] == FFRep)].copy()
			# Calculate WT_drag_contribution_summed for the subset
			WT_drag_contribution_summed = subset_df['WT_drag_contribution'].sum()
			# Calculate f_0 for the subset
			f_0 = subset_df['ref_allele_contribution'].sum()
			# Calculate d for the subset
			d = WT_drag_contribution_summed / f_0 # This is d per FF/BioRep
			effects_table.loc[(effects_table['BioRep'] == BioRep) & (effects_table['FFRep'] == FFRep), 'd'] = d
			print(subset_df)
		# Step 3: Calculate effect_size_updated for all rows -this is ei or e_new for next iteration or until the loop breaks
		effects_table['effect_size_updated'] = (2 * (effects_table['d'] * effects_table['effect_size_scaled_qpcr'] - 1)) / (1 + effects_table['adjusted_freq'])
		loops += 1
		print(loops)

	effects_table['effect_size_updated'] = effects_table['effect_size_updated'] + 1 # revert to where 1 = no effect
	effects_table = effects_table.rename(columns={'effect_size_updated': 'effect_size'}) # rename this as effect_size to be reported
	effects_table.loc[effects_table["effect_size"] < 0, "effect_size"] = 0  # update rows where effect_size < 0 to be 0
 

	# drop columns that are not needed - initial effect size (no adjustment) and effect_size_scaled_qpcr (with old adjustments) can be found in the V1 allelic effect size output file
	columns_to_drop = [
	    'initial_effect_size',
	    'effect_size_scaled_qpcr',
	    'V1_effect_size',
	    'effect_size_eold',
	    'WT_drag_contribution',
	    'ref_allele_contribution',
	    'd'
	]
	effects_table = effects_table.drop(columns=columns_to_drop)


	effects_table.to_csv(output_file, sep='\t', index=False)
	return output_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process input file for adjusting effects.')
    parser.add_argument('effects_table', type=str, help='Path to input file (TSV format)')
    parser.add_argument('output_file', type=str, help='Path to output file (TSV format)')


    args = parser.parse_args()
    adjust_effects_heterozygous_editing_V2(args.effects_table, args.output_file)