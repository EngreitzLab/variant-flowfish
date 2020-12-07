source /broad/software/scripts/useuse
reuse -q GridEngine8
use .python-3.5.1
source /seq/lincRNA/Ben/VENV_MIP/bin/activate

snakemake --cluster "/seq/lincRNA/Ben/LanderLab-EP-Prediction/src/BFF/snake-qsub" -j 100 -r -p -k -w 60

