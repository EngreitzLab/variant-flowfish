# Compute one sample two-sided T-test
# Estimate power of each variant
# Correct for multiple testing (using Benjamini-Hochberg) 


rule compute_significance_and_power:
    input:
        flat='results/summary/AllelicEffects.{replicateDirectory}.{ExperimentID}.flat.tsv.gz'
    output:
        stats='results/summary/AllelicEffectsStats.{replicateDirectory}.{ExperimentID}.tsv'
        effectCor='results/summary/effectCorrelations.pdf'
    params:
        codedir=config['codedir']
    shell:
        """
        bash -c '
                . $HOME/.bashrc 
                conda activate vff_R
                Rscript {params.codedir}/workflow/scripts/t.test.R --allelicEffectFile {input.flat} --outfile {output.stats} --repPlots {output.effectCor} --effectColumn effect_size
        """