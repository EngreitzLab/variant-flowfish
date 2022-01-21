from scripts.plot_reference_mismatches import *

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

rule plot_reference_mismatches:
    input: 
        samples = lambda wildcards:
            samplesheet.loc[samplesheet['AmpliconID']==wildcards.AmpliconID]['referenceAlleleFile']
            
    output:
        reference_plots = "results/summary/{AmpliconID}.reference_plots.pdf"
    run:
        plot_reference_mismatches(input.samples, wildcards.AmpliconID, output.reference_plots)
        
