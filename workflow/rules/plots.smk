## Plotting

rule plot_genotyping_stats:
    input:
        variantCounts="results/summary/VariantCounts.flat.tsv.gz",
        samplesheet="SampleList.snakemake.tsv"
    ## To do:  Adjust this script to use the correlation files now created in GetPCRReplicateCorrelation.R
    output: 
        "results/summary/DesiredVariants.RData"
    params:
        codedir=config['codedir']
    shell:
        "Rscript {params.codedir}/workflow/scripts/PlotVariantCounts.R --variantCounts {input.variantCounts} --samplesheet {input.samplesheet} --outbase results/summary/DesiredVariants"
