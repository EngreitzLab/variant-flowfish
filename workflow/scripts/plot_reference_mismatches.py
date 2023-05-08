import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker
import ast
import numpy as np
import seaborn as sns
import os

def plot_reference_mismatches(reference_files, amplicon, output_file):
    ref_data = []
    for file in reference_files:
        try:
            df = pd.read_csv(file, sep='\t')
            ref_data.append(df)
        except:
            # empty file
            print('%s is empty.' % file)
            continue
    ref_data = pd.concat(ref_data)

    reads_info = {}
    count = 0
    for index, row in ref_data.iterrows():
        count += 1
        if count % 10000 == 0: 
            print('processed %d/%d rows' % (count, len(ref_data)))
        for col in ['Mismatches', 'Deletions', 'Insertions']:
            for val in ast.literal_eval(row[col]):
                if val in reads_info.keys():
                    if col in reads_info[val]:
                        reads_info[val][col] += 1
                        reads_info[val][col + '#Reads'] += row['#Reads']
                    else:
                        reads_info[val][col] = 1
                        reads_info[val][col + '#Reads'] = row['#Reads']
                else:
                    reads_info[val] = {}
                    reads_info[val][col] = 1
                    reads_info[val][col + '#Reads'] = row['#Reads']

    reads_df = pd.DataFrame(reads_info).T
    reads_df = reads_df.reindex(columns=['Mismatches', 'Mismatches#Reads', 'Deletions', 'Deletions#Reads', 'Insertions', 'Insertions#Reads']) # add any missing columns
    reads_df = reads_df.fillna(0)

    # Mismatches, insertions, deletions per read
    ref_data['#Mismatches'] = ref_data['Mismatches'].apply(lambda x: len(ast.literal_eval(x)))
    ref_data['#Insertions'] = ref_data['Insertions'].apply(lambda x: len(ast.literal_eval(x)))
    ref_data['#Deletions'] = ref_data['Deletions'].apply(lambda x: len(ast.literal_eval(x)))
    
    mismatch_grouped = ref_data.groupby('#Mismatches')['#Reads'].sum().reset_index()
    mismatch_grouped['Read_freq'] = mismatch_grouped['#Reads']/(ref_data['#Reads'].sum())

    insertion_grouped = ref_data.groupby('#Insertions')['#Reads'].sum().reset_index()
    insertion_grouped['Read_freq'] = insertion_grouped['#Reads']/(ref_data['#Reads'].sum())

    deletion_grouped = ref_data.groupby('#Deletions')['#Reads'].sum().reset_index()
    deletion_grouped['Read_freq'] = deletion_grouped['#Reads']/(ref_data['#Reads'].sum())
        
    # Mismatch frequency by position  
    reads_df['MismatchFreq'] = reads_df['Mismatches#Reads']/(ref_data['#Reads'].sum())
    reads_df['InsertionFreq'] = reads_df['Insertions#Reads']/(ref_data['#Reads'].sum())
    reads_df['DeletionFreq'] = reads_df['Deletions#Reads']/(ref_data['#Reads'].sum())

    # save plots
    # locator = matplotlib.ticker.MultipleLocator(1) # force x-axis integer step size 1
    with PdfPages(output_file) as pdf:

        fig_errors, ax_errors = plt.subplots()
        ax_errors.hist(ref_data['Errors'], range=(0, 15), bins=15, density=True)
        ax_errors.set_title('%s Probability Density Histogram of Errors' % amplicon)
        ax_errors.set_xlabel('# errors (mismatches + insertions + deletions)')
        ax_errors.set_ylabel('Probability Density')
        pdf.savefig(fig_errors)
        plt.close(fig_errors)
        
        fig_mismatch, ax_mismatch = plt.subplots()
        ax_mismatch.bar(mismatch_grouped['#Mismatches'], mismatch_grouped['Read_freq'])
        ax_mismatch.set_title('%s Aligned/Reference Mismatch Frequency' % amplicon)
        ax_mismatch.set_xlabel('#Mismatches')
        ax_mismatch.set_ylabel('Read Frequency')
        # ax_mismatch.xaxis.set_major_locator(locator)
        pdf.savefig(fig_mismatch)
        plt.close(fig_mismatch)

        fig_insertion, ax_insertion = plt.subplots()
        ax_insertion.bar(insertion_grouped['#Insertions'], insertion_grouped['Read_freq'])
        ax_insertion.set_title('%s Aligned/Reference Insertion Frequency' % amplicon)
        ax_insertion.set_xlabel('#Insertions')
        ax_insertion.set_ylabel('Read Frequency')
        # ax_insertion.xaxis.set_major_locator(locator)
        pdf.savefig(fig_insertion)
        plt.close(fig_insertion)

        fig_deletion, ax_deletion = plt.subplots()
        ax_deletion.bar(deletion_grouped['#Deletions'], deletion_grouped['Read_freq'])
        ax_deletion.set_title('%s Aligned/Reference Deletion Frequency' % amplicon)
        ax_deletion.set_xlabel('#Deletions')
        ax_deletion.set_ylabel('Read Frequency')
        # ax_deletion.xaxis.set_major_locator(locator)
        pdf.savefig(fig_deletion)
        plt.close(fig_deletion)

        fig_mfreq, ax_mfreq = plt.subplots()
        ax_mfreq.bar(reads_df.index, reads_df['MismatchFreq'])
        ax_mfreq.set_title('%s Aligned/Reference Mismatch Frequency by Position' % amplicon)
        ax_mfreq.set_xlabel('Reference Position')
        ax_mfreq.set_ylabel('#Mismatch Reads/Total #Reads')
        # ax_mfreq.xaxis.set_major_locator(locator)
        ax_mfreq.ticklabel_format(axis='y', style='sci')
        pdf.savefig(fig_mfreq)
        plt.close(fig_mfreq)

        fig_ifreq, ax_ifreq = plt.subplots()
        ax_ifreq.bar(reads_df.index, reads_df['InsertionFreq'])
        ax_ifreq.set_title('%s Aligned/Reference Insertion Frequency by Position' % amplicon)
        ax_ifreq.set_xlabel('Reference Position')
        ax_ifreq.set_ylabel('#Insertion Reads/Total #Reads')
        # ax_ifreq.xaxis.set_major_locator(locator)
        ax_ifreq.ticklabel_format(axis='y', style='sci')
        pdf.savefig(fig_ifreq)
        plt.close(fig_ifreq)

        fig_dfreq, ax_dfreq = plt.subplots()
        ax_dfreq.bar(reads_df.index, reads_df['DeletionFreq'])
        ax_dfreq.set_title('%s Aligned/Reference Deletion Frequency by Position' % amplicon)
        ax_dfreq.set_xlabel('Reference Position')
        ax_dfreq.set_ylabel('#Deletion Reads/Total #Reads')
        # ax_dfreq.xaxis.set_major_locator(locator)
        ax_dfreq.ticklabel_format(axis='y', style='sci')
        pdf.savefig(fig_dfreq)
        plt.close(fig_dfreq)