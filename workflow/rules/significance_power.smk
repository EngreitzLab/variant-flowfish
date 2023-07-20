# Compute one sample two-sided T-test
# Estimate power of each variant
# Correct for multiple testing (using Benjamini-Hochberg) 


rule compute_significance_and_power:
    input:
        flat='results/summary/AllelicEffects.byExperimentRep.ExperimentIDReplicates.flat.tsv.gz'
    output:
        stats='results/summary/stats/AllelicEffectsStats.tsv',
        volcano='results/summary/stats/Volcano.pdf',
        var='results/summary/stats/modeling.pdf'
    params:
        codedir=config['codedir']
    shell:
        """
        bash -c '
                . $HOME/.bashrc 
                conda activate vff_R
                Rscript {params.codedir}/workflow/scripts/t.test.R --allelicEffectFile {input.flat} --outfile {output.stats} --varPlots {output.var} --volcanoPlots {output.volcano} --effectColumn effect_size'
        """