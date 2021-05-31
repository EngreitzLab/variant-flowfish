## Jesse Engreitz, Hank Jones, Ben Doughty
## May 4, 2021
## VFF analysis for 10-pegRNA experiments from Glen
## https://docs.google.com/spreadsheets/d/1fbP2ogkAN__4Qcfuu8Z6umAoD9H746FkFTaid2KLGKM/edit#gid=0

PROJECT=$OAK/Projects/VariantEditing/FF/210426B001-VFF/

git clone git@github.com:EngreitzLab/variant-flowfish.git


########################################################################
## Demultiplex (move this into Snakemake pipeline at some point?)

BCLDIR=$OAK/Projects/SequencingRuns/210427_SL-HDF_1225_AHL3GYBCX3/
mkdir -p $PROJECT/fastq/


## Create SampleSheet for bcl2fastq from our special VFF sample sheet format, including looking up barcode sequences
cp variant-flowfish/resources/SampleSheet.header.csv config/SampleSheet.csv
Rscript <(echo '
library(dplyr)
library(tidyr)
samplesheet <- read.delim("config/SampleList.tsv", check.names=F, stringsAsFactors=F) %>% select(SampleID,Index1,Index2)
i5 <- read.delim("variant-flowfish/resources/2Pi5IndexedPrimersBarcodes.txt") %>% select(BarcodeName,Barcode) %>% rename(Index2=BarcodeName,Index2Sequence=Barcode)
i7 <- read.delim("variant-flowfish/resources/2Pi7IndexedPrimersBarcodes.txt") %>% select(WellPosition,Barcode) %>% rename(Index1=WellPosition,Index1Sequence=Barcode)
res <- merge(samplesheet, i5, by="Index2", all.x=TRUE)
res <- merge(res, i7, by="Index1", all.x=TRUE)
res <- res %>% mutate(SampleID=ordered(SampleID, levels=samplesheet$SampleID)) %>% arrange(SampleID)
stopifnot(!any(duplicated(res$SampleID)))
stopifnot(!any(is.na(res$Index1Sequence)))
stopifnot(!any(is.na(res$Index2Sequence)))
write.table(res %>% select(SampleID,Index1Sequence,Index2Sequence) %>% rename(Sample_ID=SampleID,Index=Index1Sequence,Index2=Index2Sequence), sep=",", quote=F, row.names=F, col.names=T, file="config/SampleSheet.csv", append=TRUE)
')

## Demultiplex
bcl2fastq \
  --runfolder-dir $BCLDIR \
  --output-dir $PROJECT/fastq/ \
  --sample-sheet config/SampleSheet.csv \
  --use-bases-mask Y*,I*,I*,Y* \
  --no-lane-splitting \
  --create-fastq-for-index-reads \
  --barcode-mismatches 0 & 


######################################################################
## Run variant-flowfish pipeline

mkdir -p $PROJECT/results/ $PROJECT/log/
snakemake \
	-s variant-flowfish/workflow/Snakefile \
	--configfile config/config.json \
	--cores 1 \
	--jobs 200 \
	--cluster "sbatch -n 1 -c 1 --mem 8G -t 12:00:00 -p owners -J VFF_{rule} -o log/{rule}_{wildcards} -e log/{rule}_{wildcards}"


  mkdir -p $PROJECT/results/ $PROJECT/log/
snakemake \
  -s variant-flowfish/workflow/Snakefile \
  --configfile config/config_tmp.json \
  --cores 1 \
  --jobs 200 \
  --cluster "sbatch -n 1 -c 1 --mem 8G -t 4:00:00 -p owners -J VFF_{rule} -o log/{rule}_{wildcards} -e log/{rule}_{wildcards}"

