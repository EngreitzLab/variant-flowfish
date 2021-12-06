# Snakemake workflow: Variant-FlowFISH

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.5.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/{{cookiecutter.repo_name}}.svg?branch=master)](https://travis-ci.org/snakemake-workflows/{{cookiecutter.repo_name}})

This snakemake workflow is for analysis of Variant-FlowFISH data.


## Authors

* Ben Doughty (@bdoughty)
* Jesse Engreitz (@engreitz)
* Hank Jones
* Katherine Guo

## Description

To do

## Usage

### Step 1: Clone this github repository

[Clone](https://help.github.com/en/articles/cloning-a-repository) this to your local system, into the place where you want to perform the data analysis.

### Step 2: Install conda environment

Install conda environments using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda env create --file envs/EngreitzLab.yml
    conda env create --file envs/crispresso2_v2.2.6.yml

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 3: Set up the Sample Sheet

(Updated 11/12/21 JME)

The Sample Sheet lists all of the sequencing libraries that will be included in the analysis, and describes their relationships and groupings.
Example: https://drive.google.com/file/d/15dn5mh1MdgDYSD-yzLuXAkrItvfvt0k9/view?usp=sharing

Required columns:
    
    SampleID          Unique name for each amplicon library. (e.g., BATCH-CellLine-Sample-FFRep-PCRRep-Bin)
    AmpliconID        Name of genomic amplicon contained in the library - must match corresponding AmpliconID column in the Amplicon Table (see below)
                        Currently, this parameter (together with the Amplicon Table) controls which variants the pipeline quantifies for each sample
    Bin               Name of a FACS-sorted bin (e.g.: A B C D E F). 'All' for FlowFISH-input edited samples. 'Neg' or blank if not applicable
    PCRRep            PCR replicate number or name
    ControlForAmplicon TRUE or FALSE. Set to TRUE for unedited samples that will be used to evaluate background sequencing/PCR error rate
    EditFromGuide     [currently required, but soon not to be]: Distance of edit from guide spacer, used to control where CRISPResso looks for edits
    VFFSpikeIn        [required, but not currently used] Integer from 0 to 100 representing the percentage of unedited cells spiked into this sample

    [Experiment Keys] Provide any number of additional columns (e.g., CellLine, Guides, TestProbe) that distinguish different samples.
                        Key columns are defined as such by the 'experiment_keycols' parameter in the config file.
                        These columns will be combined to form a unique experiment key.
                        Replicates for a given unique experiment key will be combined.
                       
    [Replicate Keys]  Provide any number of additional columns (e.g., FlowFISHRep) that distinguish different experimental replicates (not including PCR replicates)
                        Replicate columns are defined as such by the 'replicate_keycols' parameter in the config file.
                        These columns will be combined to form a unique replicate id.
                        PCR replicate counts for each unique replicate key will be summed at the level of this replicate ID.
                        MLE estimates and VFF spike-in calculations will also be performed at the level of this replicate ID, 
                        then compared according to grouping of the experiment key.
                        
Optional columns:
                        
    fastqR1           If provided in the Sample Sheet, overwrites the default value (config['fastqdir']/{SampleID}_*_R1_*fastq.gz)
    fastqR2           If provided in the Sample Sheet, overwrites the default value (config['fastqdir']/{SampleID}_*_R2_*fastq.gz)
    
    Batch             Batch ID used to identify the appropriate FACS sort params file (config['sortparamsdir']/{Batch}_{SampleNumber}.csv)
    SampleNumber      FlowFISH sample number - used to identify the appropriate FACS sort params file (config['sortparamsdir']/{Batch}_{SampleNumber}.csv)
    sortParamsFile    If provided in the Sample Sheet, overwrites the default value (config['sortparamsdir']/{Batch}_{SampleNumber}.csv)
    
    Parameters that overwrite values in the Amplicon Table below:
    AmpliconSeq       Sequence of the genomic amplicon to align to. If provided in the Sample Sheet, overwrites value in the Amplicon Table for this sample.
    GuideSpacer       Spacer of the gRNA used. If provided in the Sample Sheet, overwrites value in the Amplicon Table for this sample.
    
    Parameters that overwrite values in the Variant Table below:
    VariantID         Unique readable name of the variant / allele
    MappingSequence   Genomic sequence of this variant / allele + genomic context; must match the output of CRISPResso2.    
    RefAllele         TRUE/FALSE if this is (one of) the reference alleles.  Used for plotting purposes

    Other columns can be present but are ignored.


### Step 4: Set up the Amplicon Table

(Updated 12/6/21 JME)

The Amplicon Table lists details for the genomic PCR amplicons used in the experiment.  It is optional; the information could be instead provided in the Sample Sheet.

Information from the Amplicon Table is pulled into the Sample Sheet by the 'AmpliconID' column.

Note: The AmpliconSeq should be completely spanned by the reads in the experiment — otherwise CRISPResso2 will fail.  If the sequencing reads do not entirely cover the amplicon, then adjust AmpliconSeq to match the reads.

Required columns:

    AmpliconID                  Arbitrary name of the amplicon
    AmpliconSeq                 Full genomic sequence to align to
    GuideSpacer                 Spacer sequence of the gRNA around which to quantify edits (no PAM) [not currently used]
    QuantificationWindowStart   Zero-based coordinate [) with respect to the start of the amplicon for quantifying reference allele
    QuantificationWindowEnd     Zero-based coordinate [) with respect to the start of the amplicon for quantifying reference allele 

To do:  Add additional parameters here to control the crispresso2 quantification window.


### Step 5: Provide sort params files

This file lists statistics and values derived from the FACS sort for each sample. The file names need to be named {Batch}_{SampleNumber}.csv and located in the config['sortParams'] directory, or alternatively the filename listed explicitly in the Sample Sheet in a column called 'sortParamsFile'

Required columns:

    Bin                Name of the sorted bin, needs to match "Bin" column in the Sample Sheet
    Count              Number of cells sorted into this bin
    Mean               Mean fluorescence values of cells sorted into this bin     
    Min                Minimum fluorescence value sorted into this bin (e.g., edge of the gate)
    Max                Maximum fluorescence value sorted into this bin (e.g., edge of the gate)


### Step 6: Set up the Variant Table (optional)

The Variant Table lists details for specific variants/alleles in the experiment.  It is optional, to create condensed
tables in which certain variants + alleles are named for plotting and downstream analysis.

Required columns:

    AmpliconID        Arbitrary name of the amplicon
    VariantID         Unique readable name of the variant / allele
    MappingSequence   Genomic sequence of this variant / allele + genomic context; must match the output of CRISPResso2.    
    RefAllele         TRUE/FALSE if this is (one of) the reference alleles.  Used for plotting purposes


### Step 6: Configure workflow

Edit `workflow/config.json` to point to the right files and define certain variables.

Options to control behavior of the workflow at a high level:

    genotyping_only     Set "true" to mark the workflow to process up to the point of aligning and quantifying variants, without attempting to calculate effect sizes via the maximum likelihood estimator
    single_end          Set "true" if the data is single-end as opposed to paired-end reads
    replicate_keycols   Comma-separated list of columns in the sample sheet used mark replicates (see Replicate Keys above)
    experiment_keycols  Comma-separated list of columns in the sample sheet used to mark different experiments (see Experiment Keys above)

File paths (can specify relative or absolute paths):

    sample_sheet        Path to Sample Sheet
    amplicon_info       Path to Amplicon Info table
    variant_info        Path to Variant Info table
    fastqdir            Path to directory containing FASTQ files
    sortparamsdir       Path to directory containing FACS sort parameters files [not required for genotyping_only="true"]
    codedir             Path to the variant-flowfish code (e.g. "variant-flowfish/")

Other parameters:

    crispresso_min_average_read_quality     Parameter passed to CRISPResso2 regarding minimum read quality score (e.g., 20)
    crispresso_min_single_bp_quality        Parameter passed to CRISPResso2 regarding minimum single bp quality score (e.g., 0)
    
    
### Step 7: Execute workflow

Activate the conda environment:

    conda activate EngreitzLab 

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
