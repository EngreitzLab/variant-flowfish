## Rscript to convert indices from plate format into barcode strings for demultiplexing

suppressPackageStartupMessages(library("optparse"))

option.list <- list(
  make_option("--sampleList", type="character", help="Sample sheet with columns SampleID,Index1,Index2"),
  make_option("--i5", type="character", default="variant-flowfish/resources/2Pi5IndexedPrimersBarcodes.txt", help="File with map of i5 (Index1) labels to barcodes"),
  make_option("--i7", type="character", default="variant-flowfish/resources/2Pi7IndexedPrimersBarcodes.txt", help="File with map of i7 (Index2) labels to barcodes"),
  make_option("--outfile", type="character", help="Output file"),
  make_option("--append", type="logical", default=FALSE, help="Set TRUE to append output to --outfile instead of overwriting")
  )
opt <- parse_args(OptionParser(option_list=option.list))
dput(opt)


suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

samplesheet <- read.delim(opt$sampleList, check.names=F, stringsAsFactors=F) %>% select(SampleID,Index1,Index2)
i5 <- read.delim(opt$i5)
i7 <- read.delim(opt$i7)
res <- merge(samplesheet, i5, by="Index2", all.x=TRUE)
res <- merge(res, i7, by="Index1", all.x=TRUE)
res <- res %>% mutate(SampleID=ordered(SampleID, levels=samplesheet$SampleID)) %>% arrange(SampleID)
stopifnot(!any(duplicated(res$SampleID)))
stopifnot(!any(is.na(res$Index1Sequence)))
stopifnot(!any(is.na(res$Index2Sequence)))
csv <- res %>% select(SampleID,Index1Sequence,Index2Sequence) %>% rename(Sample_ID=SampleID,Index=Index1Sequence,Index2=Index2Sequence)
write.table(csv, sep=",", quote=F, row.names=F, col.names=T, file=opt$outfile, append=opt$append)

