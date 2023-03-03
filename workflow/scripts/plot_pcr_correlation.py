import pandas as pd
from os.path import join
# from os import listdir
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy import stats
# import seaborn as sns
# sns.set(rc={'figure.figsize':(5, 5)})

def plot_pcr_correlation(pcr_replicates, output_file):
    # folder = 'results/byPCRRep'
    # p = PdfPages(join(folder, output_file))
    with PdfPages(output_file) as pdf:
        for rep in pcr_replicates:
            if rep.endswith('1.effects_vs_ref_ignoreInputBin.txt'):
                plt.clf()
                file1 = rep 
                file2 = rep[:-len('1.effects_vs_ref_ignoreInputBin.txt')] + '2.effects_vs_ref_ignoreInputBin.txt'
                rep1 = pd.read_table(file1)
                rep2 = pd.read_table(file2)
                rep1 = rep1[~rep1['RefAllele'] == True]
                rep2 = rep2[~rep2['RefAllele'] == True]
                freq_df = pd.DataFrame([np.log10(rep1['freq']), np.log10(rep2['freq'])]).T
                freq_df.columns = ['PCRRep1', 'PCRRep2']
                with pd.option_context('mode.use_inf_as_na', True):
                    freq_df.dropna(inplace=True)
                plt.scatter(freq_df['PCRRep1'], freq_df['PCRRep2'])
                plt.annotate("Pearson r = {:.3f}".format(stats.pearsonr(freq_df['PCRRep1'], freq_df['PCRRep2'])[0]), (-4, -1))
                plt.xlabel('PCR Replicate 1 (log10)')
                plt.ylabel('PCR Replicate 2 (log10)')
                plt.title(rep[len('/results/byPCRRep'):-len('-1.effects_vs_ref_ignoreInputBin.txt')] + " PCR Replicate Correlation")
                plt.xlim([-5, 0])
                plt.ylim([-5, 0])
                plt.savefig(pdf, format='pdf') 