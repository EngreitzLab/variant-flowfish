#! /usr/bin/env python

import pandas as pd
import numpy as np
import argparse
from matplotlib import pyplot as plt
plt.switch_backend('agg')

parser = argparse.ArgumentParser()
parser.add_argument("-i", dest="in_files", type=str,  nargs='+')
parser.add_argument("-o", dest="out_file", type=str)

args = parser.parse_args()

individual_effects = []

for in_file in args.in_files:
    # df = pd.read_table(in_file, usecols=['target', 'OligoID', 'mleAvg', 'sum1']).set_index('OligoID')
    df = pd.read_table(in_file, usecols=['target_no_spike', 'OligoID', 'effect_size', 'sum1_no_spike']).set_index('OligoID')
    
    if 'HEK' in in_file:
        df['class'] = 'HEK4'
    elif 'sg646' in in_file:
        df['class'] = 'sg646'
    else:
        df['class'] = 'Target'
    individual_effects.append(df)

pd.concat(individual_effects).to_csv(args.out_file, sep='\t')
