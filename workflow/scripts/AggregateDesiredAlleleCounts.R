# Jesse Engreitz
# 5/30/21
# Rscript to aggregate allele frequency information across CRISPResso runs into a table for plotting and analysis


suppressPackageStartupMessages(library("optparse"))

option.list <- list(
  make_option("--variantCounts", type="character", help="Variant count flat file (gzipped) from the variant-flowfish snakemake pipeline"),
  make_option("--variantInfo", type="character", help="File containing desired variants. Tab-delimited file containing columns AmpliconID, GuideSpacer, MappingSequence, RefAllele;  where MappingSequence matches the Aligned_Sequence output column in the Alleles frequency table in CRISPResso"),
  make_option("--minFrequency", type="numeric", default=0, help="Minimum threshold on alleles to include in merged table"),
  make_option("--outbase", type="character", default="./AlleleTable", help="Output filebase of allele frequencies (alleles x samples), which will be output with one file per AmpliconName-GuideSpacer combo")
  )
opt <- parse_args(OptionParser(option_list=option.list))
dput(opt)

if (FALSE) {
  opt <- list()
  opt$variantInfo = "../../config/AlleleList.txt"
  opt$outbase = "VariantCounts.DesiredVariants"
  opt$minFrequency = 0
  opt$variantCounts = "VariantCounts.flat.tsv.gz"
}

if (is.null(opt$outbase))
  stop("AggregateAlleleCounts: --outbase should be specified.\n")

if (is.null(opt$variantCounts) | !file.exists(opt$variantCounts)) 
  stop("AggregateAlleleCounts: --variantCounts file not found.\n")

if (is.null(opt$variantInfo) | !file.exists(opt$variantInfo))
    stop("AggregateAlleleCounts: --variantInfo file not found.\n")

if (opt$minFrequency < 0 | opt$minFrequency > 1)
  stop("AggregateAlleleCounts: --minFrequency must be between 0 and 1.\n")


suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))

# Function to compare list of target sequences to match against a longer sequence. 
# Return first match or 0 for no match.
matchSeq <- function(seq, seqList) {
  match = regmatches(seq, regexpr(paste(seqList, collapse="|"), seq)) %>% unlist()
  # use gregexpr for all matches 
  if (length(match) == 0) {
    return("None") # HELP: if something doesn't match I think the sapply I use will break 
  }
  return(match)
}


getAlleleTable <- function(countsFlat, variantInfo, minFreq) {
  
  # not sure if this is necessary, but does a similar function to the merge that was used previously
  countsFiltered <- filter(countsFlat, grepl(paste(variantInfo$MappingSequence, collapse="|"), MappingSequence)) 
  if (nrow(countsFiltered)==0) return(NULL)
  
  # apply matchSeq function to MappingSequence column and save matches in MatchSequence column
  countsFiltered["MatchSequence"] <- sapply(countsFiltered["MappingSequence"], matchSeq, seqList=variantInfo$MappingSequence)
  
  desiredCountsFlat <- countsFiltered %>% 
    select(MatchSequence,`#Reads`,`%Reads`,SampleID) %>%         
    group_by(MatchSequence,SampleID) %>% summarise(`%Reads`=sum(`%Reads`), `#Reads`=sum(`#Reads`)) %>% ungroup() %>% ## in some cases, CRISPResso reports the same Amplicon_Sequence twice ... if it aligns differently to reference
    as.data.frame()

  desiredCountsTable <- desiredCountsFlat %>% select(-`#Reads`) %>% spread("SampleID","%Reads",fill=0) %>% as.data.frame()

  if (nrow(desiredCountsTable) > 1 & ncol(desiredCountsTable) > 2) {
    ## Filter based on minimum allele frequency

    freqs <- apply(desiredCountsTable[,-1], 1, mean)
    desiredCountsTable <- desiredCountsTable[order(freqs,decreasing=T),]

    maxVal <- apply(desiredCountsTable[,-1], 1, max)
    desiredCountsTable <- desiredCountsTable[maxVal >= minFreq,]
  }

  desiredCountsFlat <- merge(variantInfo %>% mutate(origOrder=1:n()), desiredCountsFlat) %>% arrange(origOrder) %>% select(-origOrder)
  desiredCountsTable <- merge(variantInfo %>% mutate(origOrder=1:n()), desiredCountsTable) %>% arrange(origOrder) %>% select(-origOrder)

  return(list(flat=desiredCountsFlat, table=desiredCountsTable))
}


countsFlat <- read.delim(gzfile(opt$variantCounts), check.names=F, stringsAsFactors=T) %>%
            rename(MappingSequence=Aligned_Sequence) %>%
            as.data.frame()
variantInfo <- read.delim(opt$variantInfo, check.names=F, stringsAsFactors=F)

desiredCounts <- getAlleleTable(countsFlat, variantInfo, minFreq=opt$minFreq)

write.table(desiredCounts$flat, file=paste0(opt$outbase,".flat.tsv"), sep='\t', quote=F, col.names=T, row.names=F)
write.table(desiredCounts$table, file=paste0(opt$outbase,".matrix.tsv"), sep='\t', quote=F, col.names=T, row.names=F)

cat("AggregateAlleleCounts complete.\n")