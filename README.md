# Snakemake workflow: Variant-FlowFISH

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.5.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/{{cookiecutter.repo_name}}.svg?branch=master)](https://travis-ci.org/snakemake-workflows/{{cookiecutter.repo_name}})

This snakemake workflow is for analysis of Variant-FlowFISH data.


## Authors

* Ben Doughty (@bdoughty)
* Jesse Engreitz (@engreitz)

## Description

To do

## Usage

### Step 1: Clone this github repository

[Clone](https://help.github.com/en/articles/cloning-a-repository) this to your local system, into the place where you want to perform the data analysis.

### Step 2: Install conda environment

Install Snakemake and conda environment using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda env create --file envs/CRISPResso.yaml  [TODO: replace with crispresso2_env]

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 3: Set up the Sample Sheet

(Updated 5/7/21)

The Sample Sheet lists all of the sequencing libraries that will be included in the analysis, and describes their relationships and groupings.

Required columns:
    
    SampleID          Unique name for each amplicon library.
    AmpliconID        Name of genomic amplicon contained in the library - must match corresponding AmpliconID column in the Amplicon Table (see below)
    Bin               Name of a FACS-sorted bin (e.g.: A B C D E F). 'All' for input samples. 'Neg' or blank if not applicable
    PCRRep            PCR replicate number or name
    VFFSpikeIn        Integer from 0 to 100 representing the percentage of unedited cells spiked into this sample

Optional columns:

    AmpliconSeq       Sequence of the genomic amplicon to align to. If provided in the Sample Sheet, overwrites value in the Amplicon Table for this sample.
    GuideSpacer       Spacer of the gRNA used. If provided in the Sample Sheet, overwrites value in the Amplicon Table for this sample.
    [Experiment Keys] Provide any number of additional columns (e.g., CellLine) that distinguish different samples.
                        Key columns are defined as such by the 'experiment_keycols' parameter in the config file.
                        These columns will be combined to form a unique experiment key.
                        Replicates for a given unique experiment key will be combined.
    [Replicate Keys]  Provide any number of additional columns (e.g., FlowFISHRep) that distinguish different experimental replicates (not including PCR replicates)
                        Replicate columns are defined as such by the 'replicate_keycols' parameter in the config file.
                        These columns will be combined to form a unique replicate id.
                        PCR replicate counts for each unique replicate key will be summed at the level of this replicate ID.
                        MLE estimates and VFF spike-in calculations will also be performed at the level of this replicate ID, 
                        then compared according to grouping of the experiment key.
    fastqR1           If provided in the Sample Sheet, overwrites the default value (config['fastqdir']/{SampleID}_*_R1_*fastq.gz)
    fastqR2           If provided in the Sample Sheet, overwrites the default value (config['fastqdir']/{SampleID}_*_R2_*fastq.gz)


### Step 4: Set up the Amplicon Table

(Updated 5/7/21)

The Amplicon Table lists details for the genomic PCR amplicons used in the experiment.  It is optional; the information could be instead provided in the Sample Sheet.
Information from the Amplicon Table is pulled into the Sample Sheet by the 'AmpliconID' column.

Required columns:
    AmpliconID        Arbitrary name of the amplicon
    AmpliconSeq       Full genomic sequence to align to
    GuideSpacer       Spacer sequence of the gRNA around which to quantify edits (no PAM)

To do:  Add additional parameters here to control the crispresso2 quantification window.


### Step 5: Provide sort params files

To do.



### Step 6: Set up the Variant Table (optional)

The Variant Table lists details for specific variants/alleles in the experiment.  It is optional, to create condensed
tables in which certain variants + alleles are named for plotting and downstream analysis.

Required columns:
    AmpliconID        Arbitrary name of the amplicon
    VariantID         Unique readable name of the variant / allele
    MappingSequence   Genomic sequence of this variant / allele + genomic context; must match the output of CRISPResso2.    
    RefAllele         TRUE/FALSE if this is (one of) the reference alleles.  Used for plotting purposes


### Step 6: Configure workflow

[TODO]:  For now, edit `workflow/config.json` to point to the right files


### Step 7: Execute workflow

Activate the conda environment:

    conda activate EngreitzLab 
    ## TODO: Create specific environment for thie pipeline and check in yml file to workflows/envs/

Test your configuration by performing a dry-run via

    snakemake -s variant-flowfish/workflow/Snakefile --configfile config/config.json -n

Execute the workflow locally via

    snakemake -s variant-flowfish/workflow/Snakefile --configfile config/config.json -n

using `$N` cores or run it in a cluster environment (Stanford Sherlock SLURM) via

`
snakemake \
  -s variant-flowfish/workflow/Snakefile \
  --configfile config/config.json \
  --cores 1 \
  --jobs 200 \
  --cluster "sbatch -n 1 -c 1 --mem 8G -t 4:00:00 -p owners -J VFF_{rule} -o log/{rule}_{wildcards} -e log/{rule}_{wildcards}"
`

For more about cluster configuration using snakemake, see [here](https://www.sichong.site/2020/02/25/snakemake-and-slurm-how-to-manage-workflow-with-resource-constraint-on-hpc/)
