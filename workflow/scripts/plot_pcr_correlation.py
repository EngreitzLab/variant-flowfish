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

def plot_correlations_between_reps(effect_table_1, effect_table_2, title, pdf):
    bin_list = effect_table_1.filter(regex=("^[A-Z]$")).columns.tolist()
    for b in bin_list:
        freq_df = pd.DataFrame([np.log10(effect_table_1[b]), np.log10(effect_table_2[b])]).T
        freq_df.columns = ['PCRRep1', 'PCRRep2']
        with pd.option_context('mode.use_inf_as_na', True):
            freq_df.dropna(inplace=True)
        plt.figure(figsize=(6,6))
        plt.scatter(freq_df['PCRRep1'], freq_df['PCRRep2'])        
        plt.annotate("Pearson r = {:.3f}".format(stats.pearsonr(freq_df['PCRRep1'], freq_df['PCRRep2'])[0]), (-4, -1))
        plt.xlabel('PCR Replicate 1 (log10)')
        plt.ylabel('PCR Replicate 2 (log10)')
        plt.title(title + ' Bin %s' % b)
        plt.xlim([-6, 0])
        plt.ylim([-6, 0])
        plt.savefig(pdf, format='pdf') 
        plt.clf()
        
# only plotting replicate 1 and 2 correlation currently
def plot_correlations(pcr_replicates, biorep_output_file, ffrep_output_file, pcrrep_output_file):
    files = pd.DataFrame(pcr_replicates, columns=['PCRRepFile'])
    files = files.join(files.PCRRepFile.str.split(r'-([0-9]|nan)-([0-9]|nan)-([0-9]|nan)-([0-9]|nan)', \
                                    expand=True).set_axis(['prefix', 'BioRep', 'SpikeIn', 'FFRep', 'PCRRep','suffix'], \
                                                            axis=1))
    # by bio rep
    with PdfPages(biorep_output_file) as pdf:
        for pair in files.groupby(['FFRep', 'PCRRep'])['PCRRepFile'].unique():
            pair1 = pd.read_table(pair[0])
            pair2 = pd.read_table(pair[1])
            pair1 = get_freq_from_effect_table(pair1)
            pair2 = get_freq_from_effect_table(pair2)
            pair1 = pair1[~pair1['RefAllele'] == True]
            pair2 = pair2[~pair2['RefAllele'] == True]
            title = pair[0].split('.')[0].split('/')[-1] + ', ' + pair[1].split('.')[0].split('/')[-1] + '\nCorrelation'
            plot_correlations_between_reps(pair1, pair2, title, pdf)
        
    # by FF rep
    with PdfPages(ffrep_output_file) as pdf:
        for pair in files.groupby(['BioRep', 'PCRRep'])['PCRRepFile'].unique():
            pair1 = pd.read_table(pair[0])
            pair2 = pd.read_table(pair[1])
            pair1 = get_freq_from_effect_table(pair1)
            pair2 = get_freq_from_effect_table(pair2)
            pair1 = pair1[~pair1['RefAllele'] == True]
            pair2 = pair2[~pair2['RefAllele'] == True]
            title = pair[0].split('.')[0].split('/')[-1] + ', ' + pair[1].split('.')[0].split('/')[-1] + '\nCorrelation'
            plot_correlations_between_reps(pair1, pair2, title, pdf)

    # by PCR rep
    with PdfPages(pcrrep_output_file) as pdf:
        for pair in files.groupby(['BioRep', 'FFRep'])['PCRRepFile'].unique():
            pair1 = pd.read_table(pair[0])
            pair2 = pd.read_table(pair[1])
            pair1 = get_freq_from_effect_table(pair1)
            pair2 = get_freq_from_effect_table(pair2)
            pair1 = pair1[~pair1['RefAllele'] == True]
            pair2 = pair2[~pair2['RefAllele'] == True]
            title = pair[0].split('.')[0].split('/')[-1] + ', ' + pair[1].split('.')[0].split('/')[-1] + '\nCorrelation'
            plot_correlations_between_reps(pair1, pair2, title, pdf)