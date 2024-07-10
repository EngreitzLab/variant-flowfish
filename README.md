# Snakemake workflow: Variant-FlowFISH

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.5.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/{{cookiecutter.repo_name}}.svg?branch=master)](https://travis-ci.org/snakemake-workflows/{{cookiecutter.repo_name}})

This snakemake workflow is for analysis of Variant-FlowFISH data.


## Authors

* Ben Doughty (@bdoughty)
* Jesse Engreitz (@engreitz)
* Hank Jones
* Katherine Guo
* Michael Montgomery

## Description


This pipeline is configured to analyze Variant-FlowFISH and other like experiments with the capacity to analyze hundreds of variants in a single analysis. It can also be set solely to assess genome editing rates using [CRISPResso2](https://github.com/pinellolab/CRISPResso2). 

We designed this pipeline to specifically analyze genome editing efficiencies across samples/modalities, compute effect sizes for genetic variants, generate statistics concerning technical noise introduced at various steps in the molecular biological workflow, and to provide data optimized for interpretation and transferrability. Significane scores for variants are computed using 1-sample T-tests and corrected for multiple testing.


## Usage

There are several required inputs prior to installing and executing this pipeline. For ease of use, generate a subdirectory 'config/' in the directory you are performing the data analysis. Generate and place the following documents inside it.


### Input 1: Sample Sheet

The Sample Sheet lists all of the sequencing libraries that will be included in the analysis, and describes their relationships and groupings.
Example: https://drive.google.com/file/d/15dn5mh1MdgDYSD-yzLuXAkrItvfvt0k9/view?usp=sharing

Note that for FlowFISH analyses, if you provide non-sorted "Neg" samples, you will need to assign them pseudo experimental key values (AmpliconID, BioRep, Guides, etc.) from your other samples to avoid the pipeline failing to find a separate non-existant sortParams file for them.

Required columns:
    
    SampleID          Unique name for each amplicon library. (e.g., BATCH-CellLine-Sample-FFRep-PCRRep-Bin)
    AmpliconID        Name of genomic amplicon contained in the library - must match corresponding AmpliconID column in the Amplicon Table (see below)
                        Currently, this parameter (together with the Amplicon Table) controls which variants the pipeline quantifies for each sample
    
    Batch             Batch ID used to identify the appropriate FACS sort params file (config['sortparamsdir']/{Batch}_{SampleNumber}.csv)
    SampleNumber      FlowFISH sample number - used to identify the appropriate FACS sort params file (config['sortparamsdir']/{Batch}_{SampleNumber}.csv)
    
    Bin               Name of a FACS-sorted bin (e.g.: A B C D E F). 'All' for FlowFISH-input edited samples. 'Neg' or blank if not applicable
    PCRRep            PCR replicate number or name
    ControlForAmplicon TRUE or FALSE. Set to TRUE for unedited samples that will be used to evaluate background sequencing/PCR error rate
    
    [Experiment Keys] Provide any number of additional columns (e.g., CellLine, Guides, TestProbe) that distinguish different samples.
                        Key columns are defined as such by the 'experiment_keycols' parameter in the config file.
                        These columns will be combined to form a unique experiment key.
                        Replicates for a given unique experiment key will be combined.
                       
    [Replicate Keys]  Provide any number of additional columns (e.g., FlowFISHRep) that distinguish different experimental replicates (not including PCR replicates)
                        Replicate columns are defined as such by the 'replicate_keycols' parameter in the config file.
                        These columns will be combined to form a unique replicate id.
                        PCR replicate counts for each unique replicate key will be summed at the level of this replicate ID.
                        MLE estimates will also be performed at the level of this replicate ID, 
                        then compared according to grouping of the experiment key.
                        
Optional columns:
                        
    fastqR1           If provided in the Sample Sheet, overwrites the default value (config['fastqdir']/{SampleID}_*_R1_*fastq.gz)
    fastqR2           If provided in the Sample Sheet, overwrites the default value (config['fastqdir']/{SampleID}_*_R2_*fastq.gz)
    
    sortParamsFile    If provided in the Sample Sheet, overwrites the default value (config['sortparamsdir']/{Batch}_{SampleNumber}.csv)

    Parameters that overwrite values in the Amplicon Table below:
    AmpliconSeq       Sequence of the genomic amplicon to align to. If provided in the Sample Sheet, overwrites value in the Amplicon Table for this sample.
    GuideSpacer       Spacer of the gRNA used. If provided in the Sample Sheet, overwrites value in the Amplicon Table for this sample.
    
    Parameters that overwrite values in the Variant Table below:
    VariantID         Unique readable name of the variant / allele
    MappingSequence   Genomic sequence of this variant / allele + genomic context; must match the output of CRISPResso2.    
    RefAllele         TRUE/FALSE if this is (one of) the reference alleles.  Used for plotting purposes

    Other columns can be present for storing meta information but are ignored.


### Input 2: Amplicon Table


The Amplicon Table lists details for the genomic PCR amplicons used in the experiment.  It is optional; the information could be instead provided in the Sample Sheet.

Information from the Amplicon Table is pulled into the Sample Sheet by the 'AmpliconID' column.

Note: The AmpliconSeq should be completely spanned by the reads in the experiment — otherwise CRISPResso2 will fail.  If the sequencing reads do not entirely cover the amplicon, then adjust AmpliconSeq to match the reads.

Required columns:

    AmpliconID                  Arbitrary name of the amplicon
    AmpliconSeq                 Full genomic sequence to align to
    QuantificationWindowStart   Zero-based coordinate [) with respect to the start of the amplicon for quantifying reference allele
    QuantificationWindowEnd     Zero-based coordinate [) with respect to the start of the amplicon for quantifying reference allele 
    ReferenceErrorThreshold     Integer indicating how many errors (mismatch/insertion/deletion) are tolerable when inferring the reference allele.

Note: The quantification window should span the length of the amplicon you wish to assay and interpret background PCR/sequencing error ie. the region you wish to edit.

### Input 3: Sorting parameters files

This file lists statistics and values derived from the FACS sort for each sample. The file names need to be named {Batch}_{SampleNumber}.csv and located in the config['sortParams'] directory, or alternatively the filename listed explicitly in the Sample Sheet in a column called 'sortParamsFile'

Required columns:

    Name               Gate name on the cytomoter (may differ than 'Barcode' you choose to label with)
    Barcode            Name of the sorted bin, needs to match "Bin" column in the Sample Sheet
    Count              Number of cells sorted into this bin
    Mean               Mean fluorescence values of cells sorted into this bin     
    Min                Minimum fluorescence value sorted into this bin (e.g., edge of the gate)
    Max                Maximum fluorescence value sorted into this bin (e.g., edge of the gate)

Required Data:
    
    Total Sorted Population:  The sorted bins should be subpopulations of the same larger population that encompasses 
                                all the bins you sorted. For estimation purposes, we need information on the total 
                                number of cells in this population (ie. the number of cells you captured in the sort 
                                and the cells you did not sort if applicable). Label the 'Barcode' column for this 
                                population 'Total.'


### Input 4: Variant Table

The Variant Table lists details for specific variants/alleles in the experiment. 

Required columns:

    AmpliconID         Arbitrary name of the amplicon that matches AmpliconID in the provided amplicon table
    VariantID          Unique readable name of the variant / allele
    MappingSequence    Genomic sequence of this variant / allele + genomic context; We use this sequence to find the 
                        variant of interest and quantify its frequency. We typically use the variant and 3-5nt on either 
                        side of the variant to uniquely distinguish it within a given amplicon.  
    RefAllele         TRUE/FALSE if this is (one of) the reference alleles.  Used for plotting purposes

### Input 5 (optional): Guide Counts Table
TODO: Describe guide counts table skew ratio. The format is like this:

```
43936   *
25697   pegRNAsgOptiGibson-peg393-PPIF_enhancer_original-chr10:81046426:GTTAG>AGCCA
25245   pegRNAsgOptiGibson-peg359-PPIF_enhancer_original-chr10:81046401:TGGGA>AACCC
24567   pegRNAsgOptiGibson-peg406-PPIF_enhancer_original-chr10:81046441:CACCA>TCGGT
22996   pegRNAsgOptiGibson-peg388-PPIF_enhancer_original-chr10:81046421:CGTTG>TTGGG
22884   pegRNAsgOptiGibson-peg342-PPIF_enhancer_original-chr10:81046391:AGCCA>TTGGG
22720   pegRNAsgOptiGibson-peg462-PPIF_enhancer_original-chr10:81046481:GGTTT>GCAGC
22655   pegRNAsgOptiGibson-peg351-PPIF_enhancer_original-chr10:81046401:TGGGA>GCAGC
...
```

## Workflow


### Step 1: Clone this github repository

[Clone](https://help.github.com/en/articles/cloning-a-repository) this to your local system or server where you want to perform the data analysis.

### Step 2: Install conda environment

Install the "VFFenv" conda environment using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html):

We provide yml files for both the developer version (ony pins direct dependencies and defaults to minimum versions) and the release version (all dependencies frozen) of this enviroment. We recommend building the VFFenv from the release version of the yml file.
    
    mamba env create -f envs/VFFenv_release_240702.yml
    #or
    conda env create --file envs/VFFenv_dev_240702.yml

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 3: Configure workflow

Copy `example/config.json` to your newly created `config/` folder and edit the fields to point to the right files and define certain variables.

Options to control behavior of the workflow at a high level:

    genotyping_only     Set "true" to mark the workflow to process up to the point of aligning and quantifying variants, without 
                         attempting to calculate effect sizes via the maximum likelihood estimator
    single_end          Set "true" if the data is single-end as opposed to paired-end reads
    replicate_keycols   Comma-separated list of columns in the sample sheet used mark replicates (see Replicate Keys above)
    experiment_keycols  Comma-separated list of columns in the sample sheet used to mark different experiments (see Experiment Keys above)
    pooled              Set "true" if you are analyzing a pool of variants, "false" if analyzing a single variant
    
    ff_tss_guide_kd     This value is used in conjunction with the parameter below to generate a scaling factor by which to scale the
                         data by. We assume FlowFISH probesets have non-specific binding, creating a background level of fluorescence in
                         a given cell. Perform a total KD experiment of the gene of interest and measure the effects in both qPCR and FF assays. 
                         The proportional difference between the two is the scaling factor we apply to the estimated effects. Set to 1 if unknown.
    qpcr_tss_guide_kd   See above. Set to 1 if unknown.

File paths (can specify relative or absolute paths):

    sample_sheet        Path to Sample Sheet
    amplicon_info       Path to Amplicon Info table
    variant_info        Path to Variant Info table
    fastqdir            Path to directory containing FASTQ files
    sortparamsdir       Path to directory containing FACS sort parameters files [not required for genotyping_only="true"]
    codedir             Path to the variant-flowfish code (e.g. "variant-flowfish/")
    
    
### Step 4: Execute workflow

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



## Outputs


### Structure

The results generated by this workflow are partitioned in such a way to distinguish their individual contributions to the workflow. They are all housed under the aptly named '/results/' directory:

    Summary/                        Contains summary analyses for every experiment performed. This will be the most useful folder 
                                      and contains subdirectories pertaining to specific outputs to analyze.
                                
    ByPCRRep/                       Contains variant count information, mle logs, mle outputs, and PDFs quantifying effects within a given
                                      PCR replicate.
        
    ByExperimentRepCorFilter/       Contains the same as above but the analysis is performed where PCR reps are aggregated by individual 
                                      FlowFISH samples. This data also imposes a PCR correlation filter that is hard-coded. We drop samples
                                      that don't pass this threshold correlation (r = 0.8)
    ByExperimentRep/                See above minus the correlation filter
        
    aligned/                        Bam files for every sample. Includes unaligned fastq files for troubleshooting
    crispresso/                     CRISPResso2 output files
    VariantCounts/                  Count files for both reference alleles and all the variants analyzed for every fastq.
        


### Summary Subdirectories 

Within '/summary/', contain each analyses grouped by relevancy in a specific directory. Effect size analyses here are generated from folder 'ByExperimentRepCor/':

    Editing/                This folder contains information relevant for quantifying allele/variant frequencies. These data 
                              include raw frequencies, plots for visualizing WT allele frequency, and the background error rates 
                              from PCR and sequencing. The mean error rate for a given amplicon from your negative control samples or
                              unedited samples should be input into the Amplicon Table column "ReferenceErrorThreshold" to ensure
                              the WT allele is accurately quantified (failing to do so may cause inaccuracies in effect size estimations)
        
    Correlations/           Contains variant frequency correlation analyses on all the replicates and biological replicates (currently 
                             undergoing modifications). Includes quantifications of variance in frequencies at each replicate level.
                                
    Stats/                  Contains effect sizes per variant, various plots displaying this information, correlation of effect sizes,
                              variance in effect sizes.
                            
    Sequencing/             Contains read coverage and alignment count information.

        





