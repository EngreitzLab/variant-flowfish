###############################################################
## Initial example script for running Variant-FlowFISH pipeline on Stanford Sherlock
## Jesse Engreitz
## 12/14/2020

PROJECT=$OAK/Projects/VariantEditing/201214_VFFAnalysisExample; cd $PROJECT
CODEDIR=$PROJECT/variant-flowfish/

git clone git@github.com:EngreitzLab/variant-flowfish.git

CONFIG=$PROJECT/config.json
## Edit the config.json file to add path to sample sheet, etc.





###############################################################
## Run snakemake pipeline
$PROJECT/snakemake/; cd $PROJECT/snakemake/

conda activate EngreitzLab
snakemake -s "$CODEDIR/Snakefile" -j 100 -r -p \
  --directory $PROJECT/snakemake/ \
  --cluster "sbatch -J VFF -n 1 -c 1 --export=ALL --mem {cluster.memory}G --wrap" \
  --jobscript $CODEDIR/cluster/jobscript.sherlock.sh \
  --cluster-config $CODEDIR/cluster/default.json \
  -k -w 60 --configfile $CONFIG --config codedir=$CODEDIR/src/ &
