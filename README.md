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

### Step 3: Create input regions BED file

[TODO] ... e.g. Create an input 'regions' BED file, with columns `chr   start   end     name`. The snakemake script will find and score guides in these regions, and label them with the provided region name.  If desired, also collect a list of pre-designed guideRNAs (BED file with regions, and filteredGuides.bed file with guide scores)

### Step 4: Configure workflow

[TODO]

### Step 5: Execute workflow

Activate the conda environment:

    conda activate EngreitzLab 
    ## TODO: Create specific environment  

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via

`
snakemake \
  --directory results/ \
  --configfile workflow/config.json \
  --cores 1 \
  --jobs 50 \
  --cluster "sbatch -n 1 -c 1 --mem 4G -t 4:00:00 -p engreitz -J VFF_{rule} -o logs/{rule}_{wildcards}.out.txt -e logs/{rule}_{wildcards}.out.txt"
`

For more about cluster configuration using snakemake, see [here](https://www.sichong.site/2020/02/25/snakemake-and-slurm-how-to-manage-workflow-with-resource-constraint-on-hpc/)
