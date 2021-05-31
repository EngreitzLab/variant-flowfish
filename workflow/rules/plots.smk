## Plotting

rule plot_genotyping_stats:
    input:
        variantCounts="results/summary/VariantCounts.DesiredVariants.flat.tsv",
        samplesheet="SampleList.snakemake.tsv"
    output: 
        "results/summary/DesiredVariants.RData"
    params:
        codedir=config['codedir']
    shell:
        "Rscript {params.codedir}/workflow/scripts/PlotVariantCounts.R --variantCounts {input.variantCounts} --samplesheet {input.samplesheet} --outbase results/summary/DesiredVariants"
