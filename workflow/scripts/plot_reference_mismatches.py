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
        df = pd.read_csv(file, sep='\t')
        ref_data.append(df)
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
    
    # import pdb; pdb.set_trace()
    mg = plt.figure()
    plt.bar(mismatch_grouped['#Mismatches'][1:], mismatch_grouped['Read_freq'][1:])
    plt.title('%s Aligned/Reference Mismatch Frequency' % amplicon)
    plt.xlabel('#Mismatches')
    plt.ylabel('Read Frequency')
    locator = matplotlib.ticker.MultipleLocator(1)
    plt.gca().xaxis.set_major_locator(locator)

    ig = plt.figure()    
    plt.bar(insertion_grouped['#Insertions'][1:], insertion_grouped['Read_freq'][1:])
    plt.title('%s Aligned/Reference Insertion Frequency' % amplicon)
    plt.xlabel('#Insertions')
    plt.ylabel('Read Frequency')
    plt.gca().xaxis.set_major_locator(locator)

    dg = plt.figure()
    plt.bar(deletion_grouped['#Deletions'][1:], deletion_grouped['Read_freq'][1:])
    plt.title('%s Aligned/Reference Deletion Frequency' % amplicon)
    plt.xlabel('#Deletions')
    plt.ylabel('Read Frequency')
    plt.gca().xaxis.set_major_locator(locator)

    # Mismatch frequency by position  
    reads_df['MismatchFreq'] = reads_df['Mismatches#Reads']/(ref_data['#Reads'].sum())
    reads_df['InsertionFreq'] = reads_df['Insertions#Reads']/(ref_data['#Reads'].sum())
    reads_df['DeletionFreq'] = reads_df['Deletions#Reads']/(ref_data['#Reads'].sum())

    mfreq = plt.figure()
    plt.bar(reads_df.index, reads_df['MismatchFreq'])
    plt.title('%s Aligned/Reference Mismatch Frequency' % amplicon)
    plt.xlabel('Reference Position')
    plt.ylabel('#Mismatch Reads/Total #Reads')
    plt.gca().xaxis.set_major_locator(locator)

    ifreq = plt.figure()
    plt.bar(reads_df.index, reads_df['InsertionFreq'])
    plt.title('%s Aligned/Reference Insertion Frequency' % amplicon)
    plt.xlabel('Reference Position')
    plt.ylabel('#Insertion Reads/Total #Reads')
    plt.gca().xaxis.set_major_locator(locator)

    dfreq = plt.figure()
    plt.bar(reads_df.index, reads_df['DeletionFreq'])
    plt.title('%s Aligned/Reference Deletion Frequency' % amplicon)
    plt.xlabel('Reference Position')
    plt.ylabel('#Deletion Reads/Total #Reads')
    plt.gca().xaxis.set_major_locator(locator)

    pp = PdfPages(output_file)
    pp.savefig(mg)
    pp.savefig(ig)
    pp.savefig(dg)
    pp.savefig(mfreq)
    pp.savefig(ifreq)
    pp.savefig(dfreq)
    pp.close()

    # """
    # Plot the reference insertion/deletion/mismatch stats for a given amplicon.

    # Parameters
    # ----------
    # df : pandas.DataFrame
    #     DataFrame containing the reference mismatches for a given sample.
    # outdir : str
    #     Path to the output directory.
    # sample_name : str
    #     Name of the sample.
    # """
    # # Set the figure size.
    # plt.figure(figsize=(10, 6))

    # # Plot the reference mismatches.
    # sns.lineplot(x="position", y="reference_mismatch", data=df)

    # # Set the x-axis label.
    # plt.xlabel("Position")

    # # Set the y-axis label.
    # plt.ylabel("Reference mismatch")

    # # Set the title.
    # plt.title("Reference mismatch for sample {}".format(sample_name))

    # # Save the plot.
    # plt.savefig(os.path.join(outdir, "reference_mismatch_{}.png".format(sample_name)))