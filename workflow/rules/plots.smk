from scripts.plot_reference_mismatches import *

## Plotting

# Variant frequency plots
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

# Effect size plots
rule plot_allelic_effect_sizes:
	input:
		'results/{replicateDirectory}/{ExperimentIDReplicates}.effects_vs_ref.txt'
	output:
		'results/{replicateDirectory}/{ExperimentIDReplicates}.effects_vs_ref.pdf'
	params:
		codedir=config['codedir']
	shell:
		"""
		Rscript {params.codedir}/workflow/scripts/PlotMleVariantEffects.R --mleEffects {input} --outfile {output} 
		"""

rule plot_aggregate_allelic_effects:
	input:
		flat='results/summary/AllelicEffects.{replicateDirectory}.{ExperimentID}.flat.tsv.gz',
		samplesheet="SampleList.snakemake.tsv"
	output:
		'results/summary/AllelicEffects.{replicateDirectory}.{ExperimentID}.pdf'
	params:
		codedir=config['codedir']
	shell:
		"""
		Rscript {params.codedir}/workflow/scripts/PlotMleVariantEffectsAggregated.R --allelicEffectFile {input.flat} --outfile {output} --samplesheet {input.samplesheet} --effectColumn effect_size
		"""

rule plot_allelic_effect_sizes_ignoreInputBin:
	input:
		'results/{replicateDirectory}/{ExperimentIDReplicates}.effects_vs_ref_ignoreInputBin.txt'
	output:
		'results/{replicateDirectory}/{ExperimentIDReplicates}.effects_vs_ref_ignoreInputBin.pdf'
	params:
		codedir=config['codedir']
	shell:
		"""
		Rscript {params.codedir}/workflow/scripts/PlotMleVariantEffects.R --mleEffects {input} --outfile {output} 
		"""

# Reference error plots
rule plot_reference_mismatches:
    input: 
        samples = lambda wildcards:
            samplesheet.loc[samplesheet['AmpliconID']==wildcards.AmpliconID]['referenceAlleleFile']
            
    output:
        reference_plots = "results/summary/{AmpliconID}.reference_plots.pdf"
    run:
        plot_reference_mismatches(input.samples, wildcards.AmpliconID, output.reference_plots)
        
