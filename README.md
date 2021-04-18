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

    conda env create --file envs/CRISPResso.yaml  [TODO: replace]

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 3: Set up input files

[TODO]:  For the example, `cp -r /oak/stanford/groups/engreitz/Projects/VariantEditing/FF/210414_VFFPipelineTest/{ampliconinfo.txt,fastq,samplesheet.txt,sortParams} $NEWDIR`

### Step 4: Configure workflow

[TODO]:  For now, edit `workflow/config.json` to point to the right files

### Step 5: Execute workflow

Activate the conda environment:

    conda activate EngreitzLab 
    cd variant-flowfish/
    ## TODO: Create specific environment for thie pipeline and check in yml file to workflows/envs/

Test your configuration by performing a dry-run via

    snakemake --directory results/ --configfile workflow/config.json -n

Execute the workflow locally via

    snakemake --directory results/ --configfile workflow/config.json --cores 1 

using `$N` cores or run it in a cluster environment (Stanford Sherlock SLURM) via

`
snakemake \
  --directory results/ \
  --configfile workflow/config.json \
  --cores 1 \
  --jobs 50 \
  --cluster "sbatch -n 1 -c 1 --mem 4G -t 4:00:00 -p engreitz -J VFF_{rule} -o logs/{rule}_{wildcards} -e logs/{rule}_{wildcards}"
`

For more about cluster configuration using snakemake, see [here](https://www.sichong.site/2020/02/25/snakemake-and-slurm-how-to-manage-workflow-with-resource-constraint-on-hpc/)
