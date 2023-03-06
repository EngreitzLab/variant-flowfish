import pandas as pd
from os.path import join
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy import stats

def get_freq_from_effect_table(effect_table):
    effect_table = effect_table.set_index(effect_table.drop(effect_table.filter(regex=("^[A-Z]$")).columns,axis=1).columns.tolist())
    effect_table = effect_table.transform(lambda x: x / x.sum())
    effect_table.reset_index(inplace=True)
    return effect_table

# only plotting replicate 1 and 2 correlation currently
def plot_pcr_correlation(pcr_replicates, output_file):
    with PdfPages(output_file) as pdf:
        for rep in pcr_replicates:
            if rep.endswith('1.effects_vs_ref_ignoreInputBin.txt'):
                file1 = rep 
                file2 = rep[:-len('1.effects_vs_ref_ignoreInputBin.txt')] + '2.effects_vs_ref_ignoreInputBin.txt'
                rep1 = pd.read_table(file1)
                rep2 = pd.read_table(file2)
                rep1 = get_freq_from_effect_table(rep1)
                rep2 = get_freq_from_effect_table(rep2)
                rep1 = rep1[~rep1['RefAllele'] == True]
                rep2 = rep2[~rep2['RefAllele'] == True]
                bin_list = rep1.filter(regex=("^[A-Z]$")).columns.tolist()
                for b in bin_list:
                    freq_df = pd.DataFrame([np.log10(rep1[b]), np.log10(rep2[b])]).T
                    freq_df.columns = ['PCRRep1', 'PCRRep2']
                    with pd.option_context('mode.use_inf_as_na', True):
                        freq_df.dropna(inplace=True)
                    plt.figure(figsize=(6,6))
                    plt.scatter(freq_df['PCRRep1'], freq_df['PCRRep2'])
                    plt.annotate("Pearson r = {:.3f}".format(stats.pearsonr(freq_df['PCRRep1'], freq_df['PCRRep2'])[0]), (-4, -1))
                    plt.xlabel('PCR Replicate 1 Frequency (log10)')
                    plt.ylabel('PCR Replicate 2 Frequency (log10)')
                    plt.title(rep[len('/results/byPCRRep'):-len('-1.effects_vs_ref_ignoreInputBin.txt')] + " PCR Replicate Correlation Bin %s" % b)
                    plt.xlim([-6, 0])
                    plt.ylim([-6, 0])
                    plt.savefig(pdf, format='pdf') 
                    plt.clf()