source /broad/software/scripts/useuse
reuse -q GridEngine8
reuse -q Python-2.7
reuse -q Java-1.8
use .r-3.2.0-bioconductor-3.1
R_LIBS_SITE=

DESIGN=$1
COUNTS=$2
SORTDIR=$3
SORTPARAMS=$4
OUTPUT=$5
LOG=$6

Rscript /seq/lincRNA/Ben/LanderLab-EP-Prediction/src/BFF/get_allele_effect_sizes.R --designDocLocation $DESIGN --countsLocation $COUNTS --sortParamsloc ${SORTDIR}/${SORTPARAMS} --outputmle $OUTPUT --log $LOG

# Rscript /seq/lincRNA/Ben/scripts/EstimateEffectsFromBins_BFF2.R --wd $SCREEN --designDocLocation $DESIGN --countsLocation $COUNTS --repSortParamsloc ${SORTDIR}/${SORTPARAMS} --output1 $OUTPUT
