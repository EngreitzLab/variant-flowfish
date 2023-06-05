import pandas as pd
from os.path import join
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy import stats
from itertools import combinations

def get_freq_from_effect_table(effect_table):
    # making an assumption that bins are named a capital letter
    effect_table = effect_table.set_index(effect_table.drop(effect_table.filter(regex=("^[A-Z]$")).columns,axis=1).columns.tolist())
    effect_table = effect_table.transform(lambda x: x / x.sum())
    effect_table.reset_index(inplace=True)
    return effect_table

def plot_correlations_between_reps(effect_table_1, effect_table_2, title, pdf):
    bin_list = effect_table_1.filter(regex=("^[A-Z]$")).columns.tolist()
    for b in bin_list:
        freq_df = pd.DataFrame([np.log10(effect_table_1[b].reset_index(drop=True)), np.log10(effect_table_2[b].reset_index(drop=True))]).T
        freq_df.columns = ['PCRRep1', 'PCRRep2']
        with pd.option_context('mode.use_inf_as_na', True):
            freq_df.dropna(inplace=True)
        plt.figure(figsize=(6,6))
        plt.axline((0, 0), (1, 1), linewidth=1, color='k') # x=y line
        plt.scatter(freq_df['PCRRep1'], freq_df['PCRRep2'])    
        plt.annotate("Pearson r = {:.3f}".format(stats.pearsonr(freq_df['PCRRep1'], freq_df['PCRRep2'])[0]), (-4, -1))
        plt.xlabel('Replicate 1 Variant Frequency (log10)')
        plt.ylabel('Replicate 2 Variant Frequency (log10)')
        plt.title(title + ' Bin %s' % b)
        plt.xlim([-6, 0])
        plt.ylim([-6, 0])
        plt.savefig(pdf, format='pdf') 
        plt.close()

def plot_correlations_between_reps_dense(effect_table_1, effect_table_2, title, pdf):
    bin_list = effect_table_1.filter(regex=("^[A-Z]$")).columns.tolist()
    pearson_annotation = ''
    plt.figure(figsize=(6,6))
    plt.axline((0, 0), (1, 1), linewidth=1, color='k') # x=y line
    for b in bin_list:
        freq_df = pd.DataFrame([np.log10(effect_table_1[b].reset_index(drop=True)), np.log10(effect_table_2[b].reset_index(drop=True))]).T
        freq_df.columns = ['PCRRep1', 'PCRRep2']
        with pd.option_context('mode.use_inf_as_na', True):
            freq_df.dropna(inplace=True)
        plt.scatter(freq_df['PCRRep1'], freq_df['PCRRep2'], alpha=0.5, label='Bin %s' % b)        
        pearson_annotation += 'Bin %s Pearson r = {:.3f}\n'.format(stats.pearsonr(freq_df['PCRRep1'], freq_df['PCRRep2'])[0]) % b
    plt.xlim([-6, 0])
    plt.ylim([-6, 0])
    plt.legend(loc='upper left')
    plt.xlabel('Replicate 1 Variant Frequency (log10)')
    plt.ylabel('Replicate 2 Variant Frequency (log10)')
    plt.title(title)
    plt.annotate(pearson_annotation, (-3, -5.5))
    plt.savefig(pdf, format='pdf') 
    plt.close()

def plot_correlations_between_reps_effects(effect_table_1, effect_table_2, title, pdf):
    plt.figure(figsize=(6,6))
    plt.axline((0, 0), (1, 1), linewidth=1, color='k') # x=y line
    effect_table_1 = effect_table_1[(effect_table_1['sum1'] > 1000)] 
    effect_table_2 = effect_table_2[(effect_table_2['sum1'] > 1000)] 
    freq_df = pd.DataFrame([(effect_table_1['effect_size'].reset_index(drop=True)-1)*100, (effect_table_2['effect_size'].reset_index(drop=True)-1)*100]).T
    freq_df.columns = ['PCRRep1', 'PCRRep2']
    with pd.option_context('mode.use_inf_as_na', True):
        freq_df.dropna(inplace=True)
    plt.scatter(freq_df['PCRRep1'], freq_df['PCRRep2'])        
    pearson_annotation = 'Pearson r = {:.3f}\n'.format(stats.pearsonr(freq_df['PCRRep1'], freq_df['PCRRep2'])[0])#  % b
    plt.xlim([-100, 100])
    plt.ylim([-100, 100])
    plt.xlabel('PCR Replicate 1 Effect Size (%)')
    plt.ylabel('PCR Replicate 2 Effect Size (%)')
    plt.title(title + " Effect Size\nsum1 > 1000")
    plt.annotate(pearson_annotation, (-50, 50))
    plt.savefig(pdf, format='pdf') 
    plt.close()

    plt.figure(figsize=(6,6))
    plt.axline((0, 0), (1, 1), linewidth=1, color='k') # x=y line
    effect_table_1 = effect_table_1[(effect_table_1['freq'] > 0.0001)] 
    effect_table_2 = effect_table_2[(effect_table_2['freq'] > 0.0001)] 
    freq_df = pd.DataFrame([(effect_table_1['effect_size'].reset_index(drop=True)-1)*100, (effect_table_2['effect_size'].reset_index(drop=True)-1)*100]).T
    freq_df.columns = ['PCRRep1', 'PCRRep2']
    with pd.option_context('mode.use_inf_as_na', True):
        freq_df.dropna(inplace=True)
    plt.scatter(freq_df['PCRRep1'], freq_df['PCRRep2'])       
    pearson_annotation = 'Pearson r = {:.3f}\n'.format(stats.pearsonr(freq_df['PCRRep1'], freq_df['PCRRep2'])[0])#  % b
    plt.xlim([-100, 100])
    plt.ylim([-100, 100])
    plt.xlabel('Replicate 1 Effect Size (%)')
    plt.ylabel('Replicate 2 Effect Size (%)')
    plt.title(title+ " Effect Size\nsum1 > 1000 and freq > 0.0001")
    plt.annotate(pearson_annotation, (-50, 50))
    plt.savefig(pdf, format='pdf') 
    plt.close()

    plt.figure(figsize=(6,6))
    plt.axline((0, 0), (1, 1), linewidth=1, color='k') # x=y line
    effect_table_1 = effect_table_1[(effect_table_1['freq'] > 0.001)] 
    effect_table_2 = effect_table_2[(effect_table_2['freq'] > 0.001)] 
    freq_df = pd.DataFrame([(effect_table_1['effect_size'].reset_index(drop=True)-1)*100, (effect_table_2['effect_size'].reset_index(drop=True)-1)*100]).T
    freq_df.columns = ['PCRRep1', 'PCRRep2']
    with pd.option_context('mode.use_inf_as_na', True):
        freq_df.dropna(inplace=True)
    plt.scatter(freq_df['PCRRep1'], freq_df['PCRRep2'])       
    pearson_annotation = 'Pearson r = {:.3f}\n'.format(stats.pearsonr(freq_df['PCRRep1'], freq_df['PCRRep2'])[0])#  % b
    plt.xlim([-100, 100])
    plt.ylim([-100, 100])
    plt.xlabel('Replicate 1 Effect Size (%)')
    plt.ylabel('Replicate 2 Effect Size (%)')
    plt.title(title+ " Effect Size\nsum1 > 1000 and freq > 0.001")
    plt.annotate(pearson_annotation, (-50, 50))
    plt.savefig(pdf, format='pdf') 
    plt.close()
        
# only plotting replicate 1 and 2 correlation currently
def plot_pcr_correlations(pcr_replicates, biorep_output_file, ffrep_output_file, pcrrep_output_file):
    files = pd.DataFrame(pcr_replicates, columns=['PCRRepFile'])
    file_regex = r'-([0-9]*|nan)-([0-9]*|nan)-([0-9]*|nan)|\/([0-9]*|nan)-([0-9]*|nan)-([0-9]*|nan)'
    files = files.join(files.PCRRepFile.str.split(file_regex, \
                                    expand=True).dropna(axis=1).set_axis(['prefix', 'BioRep', 'FFRep', 'PCRRep','suffix'], \
                                                            axis=1))
    # by bio rep
    with PdfPages(biorep_output_file) as pdf, PdfPages(biorep_output_file.split('.')[0] + '_condensed.pdf') as pdf2, PdfPages(biorep_output_file.split('.')[0] + '_effects.pdf') as pdf3:
        for pair in files.groupby(['FFRep', 'PCRRep'])['PCRRepFile'].unique():
            if len(pair) < 2:
                print('PCR replicate pair not complete for %s ' % pair[0])
                continue
            pair1 = pd.read_table(pair[0])
            pair2 = pd.read_table(pair[1])
            pair1 = get_freq_from_effect_table(pair1)
            pair2 = get_freq_from_effect_table(pair2)
            pair1 = pair1[pair1['RefAllele'].astype('str') != 'True']
            pair2 = pair2[pair2['RefAllele'].astype('str') != 'True']
            title = pair[0].split('.')[0].split('/')[-1] + ', ' + pair[1].split('.')[0].split('/')[-1] + '\nCorrelation'
            plot_correlations_between_reps(pair1, pair2, title, pdf)
            plot_correlations_between_reps_dense(pair1, pair2, title, pdf2)
            plot_correlations_between_reps_effects(pair1, pair2, title, pdf3)
        
    # by FF rep
    with PdfPages(ffrep_output_file) as pdf, PdfPages(ffrep_output_file.split('.')[0] + '_condensed.pdf') as pdf2, PdfPages(ffrep_output_file.split('.')[0] + '_effects.pdf') as pdf3:
        for pair in files.groupby(['BioRep', 'PCRRep'])['PCRRepFile'].unique():
            if len(pair) < 2:
                print('PCR replicate pair not complete for %s ' % pair[0])
                continue
            pair1 = pd.read_table(pair[0])
            pair2 = pd.read_table(pair[1])
            pair1 = get_freq_from_effect_table(pair1)
            pair2 = get_freq_from_effect_table(pair2)
            pair1 = pair1[pair1['RefAllele'].astype('str') != 'True']
            pair2 = pair2[pair2['RefAllele'].astype('str') != 'True']
            title = pair[0].split('.')[0].split('/')[-1] + ', ' + pair[1].split('.')[0].split('/')[-1] + '\nCorrelation'
            plot_correlations_between_reps(pair1, pair2, title, pdf)
            plot_correlations_between_reps_dense(pair1, pair2, title, pdf2)
            plot_correlations_between_reps_effects(pair1, pair2, title, pdf3)

    # by PCR rep
    with PdfPages(pcrrep_output_file) as pdf, PdfPages(pcrrep_output_file.split('.')[0] + '_condensed.pdf') as pdf2, PdfPages(pcrrep_output_file.split('.')[0] + '_effects.pdf') as pdf3:
        for pair in files.groupby(['BioRep', 'FFRep'])['PCRRepFile'].unique():
            if len(pair) < 2:
                print('PCR replicate pair not complete for %s ' % pair[0])
                continue
            pair1 = pd.read_table(pair[0])
            pair2 = pd.read_table(pair[1])
            pair1 = get_freq_from_effect_table(pair1)
            pair2 = get_freq_from_effect_table(pair2)
            pair1 = pair1[pair1['RefAllele'].astype('str') != 'True']
            pair2 = pair2[pair2['RefAllele'].astype('str') != 'True']
            title = pair[0].split('.')[0].split('/')[-1] + ', ' + pair[1].split('.')[0].split('/')[-1] + '\nCorrelation'
            plot_correlations_between_reps(pair1, pair2, title, pdf)
            plot_correlations_between_reps_dense(pair1, pair2, title, pdf2)
            plot_correlations_between_reps_effects(pair1, pair2, title, pdf3)

def plot_pcr_correlations_averaged(pcr_replicates, biorep_output_file, ffrep_output_file):
    files = pd.DataFrame(pcr_replicates, columns=['PCRRepFile'])
    file_regex = r'-([0-9]*|nan)-([0-9]*|nan)-([0-9]*|nan)|\/([0-9]*|nan)-([0-9]*|nan)-([0-9]*|nan)'
    files = files.join(files.PCRRepFile.str.split(file_regex, \
                                    expand=True).dropna(axis=1).set_axis(['prefix', 'BioRep', 'FFRep', 'PCRRep','suffix'], \
                                                            axis=1))
    # by bio rep
    with PdfPages(biorep_output_file) as pdf, PdfPages(biorep_output_file.split('.')[0] + '_condensed.pdf') as pdf2, PdfPages(biorep_output_file.split('.')[0] + '_effects.pdf') as pdf3:
        bioreps = files['BioRep'].unique()
        for pair in combinations(bioreps, 2):
            br1 = files[files['BioRep'] == pair[0]]
            br1_effects_list = []
            for pcrfile in br1['PCRRepFile']:
                br1_effects_list.append(get_freq_from_effect_table(pd.read_table(pcrfile)))
            br1_effects = pd.concat(br1_effects_list)
            br1_effects_mean = br1_effects.groupby(['AmpliconID', 'MappingSequence', 'VariantID', 'RefAllele']).mean().reset_index()
            br2 = pair[1]
            br2 = files[files['BioRep'] == pair[1]]
            br2_effects_list = []
            for pcrfile in br2['PCRRepFile']:
                br2_effects_list.append(get_freq_from_effect_table(pd.read_table(pcrfile)))
            br2_effects = pd.concat(br2_effects_list)
            br2_effects_mean = br2_effects.groupby(['AmpliconID', 'MappingSequence', 'VariantID', 'RefAllele']).mean().reset_index()
            br1_effects_mean = br1_effects_mean[br1_effects_mean['RefAllele'].astype('str') != 'True']
            br2_effects_mean = br2_effects_mean[br2_effects_mean['RefAllele'].astype('str') != 'True']
            title = "BioReplicate %s, %s\nAveraged PCR Frequency Correlation" % (pair[0], pair[1])# pair[0].split('.')[0].split('/')[-1] + ', ' + pair[1].split('.')[0].split('/')[-1] + '\nCorrelation'
            plot_correlations_between_reps(br1_effects_mean, br2_effects_mean, title, pdf)
            plot_correlations_between_reps_dense(br1_effects_mean, br2_effects_mean, title, pdf2)
            plot_correlations_between_reps_effects(br1_effects_mean, br2_effects_mean, title, pdf3)
        
        # by FF rep
        with PdfPages(ffrep_output_file) as pdf, PdfPages(ffrep_output_file.split('.')[0] + '_condensed.pdf') as pdf2, PdfPages(ffrep_output_file.split('.')[0] + '_effects.pdf') as pdf3:
            for b in bioreps:
                ff_reps = files[files['BioRep'] == b]['FFRep'].unique()
                for pair in combinations(ff_reps, 2):
                    ff1 = files[(files['BioRep'] == b) & (files['FFRep'] == pair[0])]
                    ff1_effects_list = []
                    for pcrfile in ff1['PCRRepFile']:
                        ff1_effects_list.append(get_freq_from_effect_table(pd.read_table(pcrfile)))
                    ff1_effects = pd.concat(ff1_effects_list)
                    ff1_effects_mean = ff1_effects.groupby(['AmpliconID', 'MappingSequence', 'VariantID', 'RefAllele']).mean().reset_index()
                    ff1_effects_mean = ff1_effects_mean[ff1_effects_mean['RefAllele'].astype('str') != 'True']
                    ff2 = files[(files['BioRep'] == b) & (files['FFRep'] == pair[1])]
                    ff2_effects_list = []
                    for pcrfile in ff2['PCRRepFile']:
                        ff2_effects_list.append(get_freq_from_effect_table(pd.read_table(pcrfile)))
                    ff2_effects = pd.concat(ff2_effects_list)
                    ff2_effects_mean = ff2_effects.groupby(['AmpliconID', 'MappingSequence', 'VariantID', 'RefAllele']).mean().reset_index()
                    ff2_effects_mean = ff2_effects_mean[ff2_effects_mean['RefAllele'].astype('str') != 'True']
                    title = "BioReplicate %s, FFReplicate %s, %s\nAveraged PCR Frequency Correlation" % (b, pair[0], pair[1])
                    plot_correlations_between_reps(ff1_effects_mean, ff2_effects_mean, title, pdf)
                    plot_correlations_between_reps_dense(ff1_effects_mean, ff2_effects_mean, title, pdf2)
                    plot_correlations_between_reps_effects(ff1_effects_mean, ff2_effects_mean, title, pdf3)
