#! /usr/bin/env python

import pandas as pd
import numpy as np
import argparse
from matplotlib import pyplot as plt
plt.switch_backend('agg')
import seaborn as sns
from pathlib import Path


parser = argparse.ArgumentParser()

parser.add_argument("-i", dest="in_files", type=str,  nargs='+')
parser.add_argument("-s", dest="scatter", type=str)
parser.add_argument("-b", dest="boxplot", type=str)
parser.add_argument("-t", dest="hist", type=str)
parser.add_argument("-a", dest="avg", type=str)

args = parser.parse_args()

screen_effects = [pd.read_table(f) for f in args.in_files]

assert len(screen_effects) == 2

rep1, rep2 = screen_effects[0], screen_effects[1]

merged = rep1.merge(rep2, on='OligoID', suffixes=['_1', '_2'], how='inner').set_index('OligoID')
merged['counts_1'] = merged['sum1_1'].rank(pct=True)
merged['counts_2'] = merged['sum1_2'].rank(pct=True)

# threshold
filtered = merged.loc[(merged.counts_1 > 0.01) & (merged.counts_2 > 0.01)].copy()

# scatter pairs of effects
fg = sns.FacetGrid(data=filtered, hue='class_1', aspect=1.61)
fg.map(plt.scatter, 'mleAvg_1', 'mleAvg_2').add_legend()
plt.savefig(args.scatter)
plt.clf()

# average effects
filtered['mleAvg'] = filtered.filter(like='mleAvg').mean(axis=1)

# hist average
plt.hist(filtered.loc[filtered['class_1'] != 'Target', 'mleAvg'], bins='auto', label='Controls')
mles = filtered.loc[(filtered.class_1 == 'Target') & (filtered.target_1 != 'negative_control'), 'mleAvg'].values
plt.scatter(mles, np.zeros(mles.shape), c='r', label='Target')
plt.xlabel('mleAvg')
plt.ylabel('Frequency')
plt.legend()
plt.savefig(args.hist)
plt.clf()

# boxplot average
bp = filtered.boxplot(column='mleAvg', by='class_1', grid=False)

categories = ['HEK4', 'Target', 'sg646']

for i in [0, 1, 2]:
    y = filtered.mleAvg[filtered['class_1']==categories[i]].dropna()
    x = np.random.normal(i, 0.025, size=len(y))
    plt.plot(x+1, y, 'r.', alpha=0.4)
plt.title('MLE Distribution across classes')
plt.ylabel('MLE')
plt.suptitle('')
plt.savefig(args.boxplot)


# write average
filtered.to_csv(args.avg, sep='\t')

