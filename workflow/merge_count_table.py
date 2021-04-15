#! /usr/bin/env python

import pandas as pd
import numpy as np
import argparse
from matplotlib import pyplot as plt
plt.switch_backend('agg')

parser = argparse.ArgumentParser()
parser.add_argument("-i", dest="in_files", type=str,  nargs='+')
parser.add_argument("-o", dest="out_file", type=str)
parser.add_argument("-p", dest='plot_file', type=str)
parser.add_argument("-s", dest="sortparam", type=str)
parser.add_argument("-l", dest="log", type=str)
parser.add_argument("-t", dest="threshold", type=int, default=200)

args = parser.parse_args()

# how many reads a column needs to be viable... how to determine?
read_threshold = args.threshold

# constants
bins = ['A', 'B', 'C', 'D', 'E', 'F']
colors = {'A': 'red', 'B': 'orange', 'C': 'yellow', 'D': 'green', 'E': 'blue', 'F': 'purple'}

def scatter_replicates(count_table):
    percentage_table = count_table.div(count_table.sum(axis=0))*100

    bounds = [100, 1.1*np.sort(percentage_table.values.reshape(-1))[-13]]

    fig, axes = plt.subplots(2, figsize=(5,10))

    for idx, ax in enumerate(axes):
        ax.plot((0, bounds[idx]), (0, bounds[idx]))

        for b in bins:
            binned = percentage_table.filter(regex='.*-{}$'.format(b))
            binned.plot.scatter(binned.columns[0], binned.columns[1], ax=ax, color=colors[b], label=b)

        ax.set(adjustable='box-forced', aspect='equal')

        ax.set_xlim([0, bounds[idx]])
        ax.set_ylim([0, bounds[idx]])
        ax.legend()
        ax.set_xlabel('Rep1')
        ax.set_ylabel('Rep2')
        ax.set_title('Allele frequency in PCR replicates')
    plt.savefig(args.plot_file)

with open(args.log, 'w') as log:


    # combine input tables
    merged = pd.read_table(args.in_files.pop(0)).set_index('MappingSequence')
    for df in args.in_files:
        merged = merged.join(pd.read_table(df).set_index('MappingSequence'), how='outer')

    merged = merged.fillna(0)

    # do plotting (before dropping columns)
    scatter_replicates(merged)

    # drop columns (bins) with too few reads
    to_drop = merged.sum(axis=0) < read_threshold
    filtered = merged.loc[:,~to_drop]

    log.write('Wells dropped\n')
    for dropped in to_drop.loc[to_drop].iteritems():
        log.write(dropped[0]+'\n')

    #filtered = merged.loc[:, (merged.sum(axis=0) > read_threshold)]

    #print(filtered.columns)

    # convert to percentages
    filtered = filtered.div(filtered.sum(axis=0))*100

    # drop rows (alleles) with too few reads
    filtered = filtered.loc[filtered.mean(axis=1) > 0.05]

    for b in bins:
        to_average = filtered.filter(regex='.*-{}$'.format(b))
        if to_average.empty:
            filtered['{}-{}'.format(args.sortparam, b)] = 0
        else:
            filtered['{}-{}'.format(args.sortparam, b)] = to_average.mean(axis=1)

    filtered[['{}-{}'.format(args.sortparam, b) for b in bins]].to_csv(args.out_file, sep='\t')

    
