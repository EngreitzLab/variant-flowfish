#! /usr/bin/env python

import pandas as pd
import numpy as np
import argparse
from matplotlib import pyplot as plt
plt.switch_backend('agg')
from scipy.optimize import root

# QUESTIONS: Why is this sometimes returning a negative value? How do we convert these numbers into effect sizes? Do we need to convert to real space differently (i.e. using SD)? 

def correct_with_spikein(no_spike_file, spike_file, corrected_file, delta_rho, log):
    # read in files
    no_spike = pd.read_table(no_spike_file).set_index('OligoID')
    spike = pd.read_table(spike_file).set_index('OligoID')
    
    # add mean column
    spike['mean'] = np.power(10, spike['logMean']+.5*spike['logSD']**2)
    no_spike['mean'] = np.power(10, no_spike['logMean']+.5*no_spike['logSD']**2)

    # merge
    merged = no_spike[['target', 'sum1', 'logMean', 'logSD', 'mean']].join(spike[['target', 'sum1', 'logMean', 'logSD', 'mean']], lsuffix='_no_spike', rsuffix='_spike', how='inner')

    # calculate scaling factor (between experiments)
    merged['alpha'] = merged['mean_spike'] / merged['mean_no_spike']

    # figure out fraction/weight to assign each allele
    merged['weights'] = merged['sum1_no_spike'] / merged.loc[merged.target_spike == '*', 'sum1_no_spike'].sum() 

    merged['freq'] = merged['sum1_no_spike'] / merged['sum1_no_spike'].sum() 

    # calculate average scaling factor
    alpha = (merged.loc[merged.target_spike == '*', 'alpha'] * merged.loc[merged.target_spike == '*', 'weights']).sum()
    log.write('alpha: {}\n'.format(alpha))
    
    # correct between experiments
    merged['mean_spike_corrected'] = merged.mean_spike / alpha

    # get WT expressions
    beta_0 = merged.loc[merged.target_no_spike == 'negative_control', 'mean_no_spike'].values[0]
    beta_0_spike = merged.loc[merged.target_no_spike == 'negative_control', 'mean_spike_corrected'].values[0]
    log.write('beta_0: {}\n'.format(beta_0))
    log.write('beta_0_spike: {}\n'.format(beta_0_spike))

    # get WT freq
    f_0 = merged.loc[merged.target_no_spike == 'negative_control', 'freq'].values[0]
    log.write('f_0: {}\n'.format(f_0))

    # get weighted average of betas
    ### DON'T USE WEIGHTS HERE!!! THEY SUM TO > 1, SINCE THEY WERE JUST FOR CALCULATING THE ALPHA
    merged['weighted_eff'] = merged['freq'] * merged['mean_no_spike'] #np.power(10, merged['logMean_no_spike'])
    beta_avg = merged['weighted_eff'].sum()
    log.write('beta_avg: {}\n'.format(beta_avg))

    # calculate deltaB
    delta_beta = beta_0 - beta_0_spike # np.power(10, merged.loc[merged.target_no_spike == 'negative_control', 'logMean_spike'].values[0])
    log.write('delta_beta: {}\n'.format(delta_beta))

    # define optimization objective and solve for rho
    to_optimize = lambda rho: delta_beta - delta_rho * (beta_avg - rho * beta_0) / ((f_0 + delta_rho) * (2 * f_0 - 3 * f_0 * rho + rho ** 2))
    rho = root(to_optimize, [.5])['x'][0]
    log.write('rho: {}\n'.format(rho))

    # solve for e_avg (without rho)
    eps_avg = delta_beta * (f_0 + delta_rho) * f_0 / delta_rho
    log.write('eps_avg: {}\n'.format(eps_avg))

    # correct remaining 
    merged['eps'] = merged['mean_no_spike'] - eps_avg
    merged.loc[merged.target_no_spike == 'negative_control', 'eps'] = merged.loc[merged.target_no_spike == 'negative_control', 'mean_no_spike'] - min(1, (1 - rho / f_0)) * eps_avg

    # normalize to WT
    merged['effect_size'] = merged['eps'] / merged.loc[merged.target_no_spike == 'negative_control', 'eps'].values[0]

    # write out
    # merged.reset_index()[['OligoID', 'target_no_spike', 'effect_size']].to_csv(args.corrected, index=None, sep='\t')                                                                                                                                                 
    # pd.read_table(args.no_spike).to_csv(args.corrected, index=None, sep='\t')
    merged['effect_size'] = merged['mean_no_spike'] / merged.loc[merged.target_no_spike == 'negative_control', 'mean_no_spike'].values[0]
    merged.reset_index()[['OligoID', 'target_no_spike', 'effect_size', 'sum1_no_spike']].to_csv(corrected_file, index=None, sep='\t')





    #deltaB = delta_rho*(beta_avg - rho * beta_0) / ((f_0 + delta_rho) * (2*f_0 - 3*f_0*rho + rho ** 2))

    # no_spike = pd.read_table(no_spike_file).set_index('OligoID')[['target', 'sum1', 'logMean', 'logSD']]
    # spike = pd.read_table(spike_file).set_index('OligoID')[['target', 'sum1', 'logMean', 'logSD']]

    # merged = no_spike.join(spike, lsuffix='_no_spike', rsuffix='_spike', how='inner')

    # beta_0 = np.power(10, merged.loc[merged.target_no_spike == 'negative_control', 'logMean_no_spike'].values[0])

    # merged['freq_no_spike'] = merged['sum1_no_spike'] / merged['sum1_no_spike'].sum()
    # f_0 = merged.loc[merged.target_no_spike == 'negative_control', 'freq_no_spike'].values[0]

    # merged['weighted_eff'] = merged['freq_no_spike'] * np.power(10, merged['logMean_no_spike'])
    # beta_avg = merged['weighted_eff'].sum()

    # delta_beta = beta_0 - np.power(10, merged.loc[merged.target_no_spike == 'negative_control', 'logMean_spike'].values[0])

    



parser = argparse.ArgumentParser()
parser.add_argument("-n", dest="no_spike", type=str)
parser.add_argument("-s", dest="spike", type=str)
parser.add_argument("-c", dest='corrected', type=str)
parser.add_argument("-d", "--delta_rho", dest="delta_rho", type=float, default=1)
parser.add_argument("-l", dest="log", type=str)

args = parser.parse_args()

with open(args.log, 'w') as log:
    correct_with_spikein(args.no_spike, args.spike, args.corrected, args.delta_rho, log)



