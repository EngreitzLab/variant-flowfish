# Jesse Engreitz
# 5/30/21


suppressPackageStartupMessages(library("optparse"))

option.list <- list(
  make_option("--variantCounts", type="character", help="Desired variant count flat file"),
  make_option("--samplesheet", type="character", help="Snakemake sample sheet including SampleID info"),
  make_option("--minCorrelation", type="numeric", default=-1, help="PCR replicates with less than this correlation will be flagged for removal"),
  make_option("--correlationFile", type="character", help="Output file for PCR replicate correlation table"),
  make_option("--lowCorSamplesFile", type="character", help="Output file to list samples flagged with low replicate correlations"),
  make_option("--cvFile", type="character", help="Output file for coefficient of variation table")
  )
opt <- parse_args(OptionParser(option_list=option.list))
dput(opt)

if (is.null(opt$correlationFile) | is.null(opt$cvFile))
  stop("GetPCRReplicateCorrelation: --correlationFile and --cvFile should be specified.\n")

if (is.null(opt$variantCounts) | !file.exists(opt$variantCounts))
  stop("GetPCRReplicateCorrelation: --variantCounts file not found.\n")

if (is.null(opt$samplesheet) | !file.exists(opt$samplesheet))
    stop("GetPCRReplicateCorrelation: --samplesheet file not found.\n")

if (opt$minCorrelation < -1 | opt$minCorrelation > 1)
  stop("GetPCRReplicateCorrelation: --minCorrelation must be [-1,1]")


suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))

countsFlat <- read.delim(opt$variantCounts, check.names=F, stringsAsFactors=F)
samplesheet <- read.delim(opt$samplesheet, check.names=F, stringsAsFactors=F)

binList <- unique(samplesheet$Bin)
binList <- binList[!(binList %in% c("All","Neg",""))]


############################################
## Calculate correlation and coefficient of variation tables


flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
    )
}


getLowCorrelationSamples <- function(cormat, minCorrelation) {
  ## Flag samples from the correlation matrix if they show low correlations with other samples
  ## Recursively remove the worst samples until either the remaining matrix is above the minimum threshold,
  ## or there is only 1 sample left in which case return the name of this sample too
  if (nrow(cormat) == 1) return(rownames(cormat))

  lowCorrelations <- apply(cormat, 1, function(row) sum(row <= minCorrelation))
  meanCorrelation <- apply(cormat, 1, mean)
  lowCorrelations[is.na(lowCorrelations)] <- 0
  meanCorrelation[is.na(meanCorrelation)] <- 0
  if (all(lowCorrelations == 0)) return(c())

  i <- which.max(lowCorrelations * (1-meanCorrelation))[1]
  worstSample <- rownames(cormat)[i]
  return(c(worstSample, getLowCorrelationSamples(cormat[-i,-i,drop=F], minCorrelation)))
}


getPCRReplicateCorrelations <- function(countsFlat, samplesheet, includeRef=FALSE, minCorrelation=-1) {
  ## Return correlations among all pairs of PCR replicates, considering editing rates of non-reference desired alleles
  results <- list()
  lowCorSamples <- list()
  samplesheet <- samplesheet %>% mutate(Grouping=paste0(ExperimentIDReplicates,Bin))
  for (group in unique(samplesheet$Grouping)) {
    currSamples <- samplesheet %>% filter(Grouping == group)

    currCounts <- countsFlat %>%
      filter(SampleID %in% currSamples$SampleID) %>%
      filter(includeRef | RefAllele == "False") %>%  # it's a string not a binary
      select(SampleID,VariantID,`%Reads`) %>%
      spread(SampleID,`%Reads`,fill=0) %>%
      select(-VariantID) %>%
      as.matrix()

    correlations <- cor(currCounts)
    if (nrow(correlations) > 1) {
      lowCorSamples[[group]] <- getLowCorrelationSamples(correlations, minCorrelation)
      cor.df <- flattenCorrMatrix(correlations)
      results[[group]] <- cor.df
    }
  }
  flat <- do.call(rbind, results)
  return(list(cor=flat, lowCorSamples=unlist(lowCorSamples)))
}


getPCRReplicateVariantCV <- function(countsFlat, samplesheet, includeRef=FALSE) {
  ## Return coefficient of variation for each variant among all pairs of PCR replicates — so that we can see the degree of variance as a function of allele frequency
  results <- list()
  samplesheet <- samplesheet %>% mutate(Grouping=paste0(ExperimentIDReplicates,Bin))
  for (group in unique(samplesheet$Grouping)) {
    currSamples <- samplesheet %>% filter(Grouping == group)
    if (nrow(currSamples) >= 2) {
      tryCatch({
        currCounts <- countsFlat %>%
          filter(SampleID %in% currSamples$SampleID) %>%
          filter(includeRef | RefAllele == "False") %>%
          select(SampleID,VariantID,`%Reads`) %>%
          spread(SampleID,`%Reads`,fill=0) %>%
          select(-VariantID) %>%
          as.matrix()

        cv <- t(apply(currCounts, 1, function(row) return(c(mean(row), sd(row), sd(row)/mean(row)*100))))
        if (ncol(cv) == 3)
          results[[group]] <- cv
      }, error = function(e) print("Failed on group = ", group))
    }
  }
  flat <- data.frame(do.call(rbind, results))
  colnames(flat) <- c("mean", "sd", "CV")
  return(flat)
}


cor.df <- getPCRReplicateCorrelations(countsFlat, samplesheet, minCorrelation=opt$minCorrelation)
vcv <- getPCRReplicateVariantCV(countsFlat, samplesheet)

write.table(cor.df$cor, file=opt$correlationFile, sep='\t', quote=F, col.names=T, row.names=F)
write.table(cor.df$lowCorSamples, file=opt$lowCorSamplesFile, sep='\t', quote=F, col.names=F, row.names=F)
write.table(vcv, file=opt$cvFile, sep='\t', quote=F, col.names=T, row.names=F)
