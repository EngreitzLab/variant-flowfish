### Demultiplex FASTQ files from the Sample Sheet
## NOT FINISHED


#BCLDIR=$OAK/Projects/SequencingRuns/210427_SL-HDF_1225_AHL3GYBCX3/
#210427_SL-HDF_1225_AHL3GYBCX3

  "bcl2fastq \
  --runfolder-dir $PROJECT/210423_NS500735_0187_AHVLHKBGXG/ \
  --output-dir $PROJECT/fastq/ \
  --sample-sheet SampleSheet.csv \
  --use-bases-mask Y*,I*,I*,Y* \
  --no-lane-splitting \
  --create-fastq-for-index-reads \
  --barcode-mismatches 0 \
  --mask-short-adapter-reads 8" 
