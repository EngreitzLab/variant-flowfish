from scripts.plot_reference_mismatches import *
from scripts.plot_pcr_correlation import *

## Plotting

# CRISPResso aggregate plot
rule plot_crispresso_aggregate_reads:
	input:
		CRISPRessoAggFolder="results/crispresso/CRISPRessoAggregate_on_Aggregate/",
	output:
		"results/summary/crispresso_aggregate_reads.pdf"
	params:
		codedir=config['codedir']
	shell:
		"""
		Rscript {params.codedir}/workflow/scripts/PlotCRISPRessoReads.R --CRISPRessoAggFolder {input.CRISPRessoAggFolder} --genotyping_only FALSE --outfile {output} 
		"""

rule plot_crispresso_aggregate_reads_genotyping_only:
	input:
		CRISPRessoAggFolder="results/crispresso/CRISPRessoAggregate_on_Aggregate/",
	output:
		"results/summary/crispresso_aggregate_reads_genotyping_only.pdf"
	params:
		codedir=config['codedir']
	shell:
		"""
		Rscript {params.codedir}/workflow/scripts/PlotCRISPRessoReads.R --CRISPRessoAggFolder {input.CRISPRessoAggFolder} --genotyping_only TRUE --outfile {output} 
		"""

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

# PCR correlation plots
if not genotyping_only:
	rule plot_PCR_correlations:
		input:
			pcr_replicates = samplesheet['ExperimentIDPCRRep_BinCounts'].str.rstrip('.bin_counts.txt').unique() + '.effects_vs_ref_ignoreInputBin.txt'
		output:
			biorep_correlation_plots = "results/summary/correlation_plots/BioReplicatePCRCorrelations.pdf",
			ffrep_correlation_plots = "results/summary/correlation_plots/FFReplicatePCRCorrelations.pdf",
			pcr_correlation_plots = "results/summary/correlation_plots/PCRReplicateCorrelations.pdf"
		run:
			plot_pcr_correlations(input.pcr_replicates, output.biorep_correlation_plots, output.ffrep_correlation_plots, output.pcr_correlation_plots)

	rule plot_experiment_correlations:
		input:
			pcr_replicates = samplesheet['ExperimentIDPCRRep_BinCounts'].str.rstrip('.bin_counts.txt').unique() + '.effects_vs_ref_ignoreInputBin.txt'
		output:
			biorep_correlation_plots = "results/summary/correlation_plots/BioReplicateCorrelations.pdf",
			ffrep_correlation_plots = "results/summary/correlation_plots/FFReplicateCorrelations.pdf",
		run:
			plot_pcr_correlations_averaged(input.pcr_replicates, output.biorep_correlation_plots, output.ffrep_correlation_plots)



# Reference error plots
rule plot_reference_mismatches:
    input: 
        samples = lambda wildcards:
            samplesheet.loc[samplesheet['AmpliconID']==wildcards.AmpliconID]['referenceAlleleFile']
            
    output:
        reference_plots = "results/summary/{AmpliconID}.reference_plots.pdf"
    run:
        plot_reference_mismatches(input.samples, wildcards.AmpliconID, output.reference_plots)



# rule plot_variant_stats:
# 	input:
# 		stats='results/summary/AllelicEffectsStats.tsv',
# 		allEffects='results/summary/AllelicEffects.byExperimentRep.ExperimentIDReplicates.flat.tsv.gz'
# 	output: 
#   		power='results/summary/powerPlots.pdf',
#   		volcano='results/summary/volcanoEffects.pdf'
# 	params:
#   		codedir=config['codedir'],
#   		reps=config['reps']
# 	shell:
# 			"""
# 			bash -c '
# 				. $HOME/.bashrc 
# 				conda activate VFFenv
# 				Rscript {params.codedir}/workflow/scripts/VolcanoPlot.R --variantStatsFile {input.stats} --allEffects {input.allEffects} --powerPlots {output.power} --volcanoPlots {output.volcano} --rep {params.reps}'
# 			"""









        
