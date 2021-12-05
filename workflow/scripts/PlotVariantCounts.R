# Jesse Engreitz
# 5/30/21
# Rscript to aggregate allele frequency information across CRISPResso runs into a table for plotting and analysis


suppressPackageStartupMessages(library("optparse"))

option.list <- list(
  make_option("--variantCounts", type="character", help="Desired variant count flat file"),
  make_option("--samplesheet", type="character", help="Snakemake sample sheet including SampleID info"),
  #make_option("--groupby", type="character", default=NULL, help="column to group"),
  #make_option("--experimentKeyCols", type="character", help="Comma-separated list of experimental key columns in the sample sheet, as specified in snakemake config file"),
  #make_option("--replicateKeyCols", type="character", help="Comma-separated list of replicate key columns in the sample sheet, as specified in snakemake config file"),
  #make_option("--variantInfo", type="character", help="File containing desired variants. Tab-delimited file containing columns AmpliconID, GuideSpacer, MappingSequence, RefAllele;  where MappingSequence matches the Aligned_Sequence output column in the Alleles frequency table in CRISPResso"),
  make_option("--outbase", type="character", default="./DesiredVariant", help="Output filebase of plot file")
  )
opt <- parse_args(OptionParser(option_list=option.list))
dput(opt)

if (FALSE) {
  ## For testing purposes
  opt <- list()
  #opt$variantInfo = "../../config/AlleleList.txt"
  opt$outbase = "results/summary/DesiredVariants"
  opt$samplesheet = "SampleList.snakemake.tsv"
  opt$variantCounts = "results/summary/VariantCounts.DesiredVariants.flat.tsv"
}


if (is.null(opt$outbase))
  stop("PlotVariantCounts: --outbase should be specified.\n")

if (is.null(opt$variantCounts) | !file.exists(opt$variantCounts))
  stop("PlotVariantCounts: --variantCounts file not found.\n")

if (is.null(opt$samplesheet) | !file.exists(opt$samplesheet))
    stop("PlotVariantCounts: --samplesheet file not found.\n")

#experimentKeys <- strsplit(opt$experimentKeyCols,",")[[1]]
#if (length(experimentKeys) == 0)
#  stop("Could not parse --experimentKeyCols")

#replicateKeys <- strsplit(opt$replicateKeyCols,",")[[1]]
#if (length(replicateKeys) == 0)
#  stop("Could not parse --replicateKeyCols")


suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))


countsFlat <- read.delim(opt$variantCounts, check.names=F, stringsAsFactors=F)
samplesheet <- read.delim(opt$samplesheet, check.names=F, stringsAsFactors=F)

binList <- unique(samplesheet$Bin)
binList <- binList[!(binList %in% c("All","Neg",""))]



############################################
## Plot overall edited rate in a stacked barplot

getStackedBarplot <- function(countsFlat, samples, group="ExperimentIDPCRRep", fill="VariantID", includeRef=FALSE, plotNReads=FALSE) {
  counts <- countsFlat %>%
            filter(SampleID %in% samples$SampleID) %>%
            filter(includeRef | RefAllele == "False") %>%
            merge(samples %>% select("SampleID", group, "ControlForAmplicon","CellLine")) %>%
            dplyr:::rename(Frequency="%Reads", nReads="#Reads") %>%
            mutate(Edited=ordered(ControlForAmplicon, levels=c(FALSE,TRUE), labels=c("Edited","Unedited"))) %>%
            as.data.frame()

  y <- ifelse(plotNReads, "nReads", "Frequency")
  ylab <- ifelse(plotNReads, "Variant Read Count (#)", "Variant Frequency (%)")
  p <- ggplot(counts, aes_string(x=group, y=y, fill=fill)) + geom_col()
  p <- p + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6)) + ylab(ylab)
  if (plotNReads) p <- p + scale_y_continuous(trans='log10')
  p <- p + theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
                 legend.title = element_text(size=9), #change legend title font size
                 legend.text = element_text(size=8)) #change legend text font size
  p <- p + facet_grid(cols=vars(Edited), scales = "free", space = "free")
  return(p)
}


pdf(file=paste0(opt$outbase, ".totalEditing.stackedBarplots.pdf"), width=9, height=8)
samplesInput <- samplesheet %>% filter(Bin == "All" | ControlForAmplicon)
for (amplicon in unique(samplesheet$AmpliconID)) {
  currSamples <- samplesheet %>% filter(AmpliconID == amplicon)
  p <- getStackedBarplot(countsFlat, currSamples) + ggtitle(paste0("AmpliconID==",amplicon))
  print(p)
  p <- getStackedBarplot(countsFlat, currSamples, plotNReads=TRUE) + ggtitle(paste0("AmpliconID==",amplicon))
  print(p)
}
invisible(dev.off())


############################################
## Plot overall reference allele rate in a side-by-side barplot

getRefStackedBarplot <- function(countsFlat, samples, group="ExperimentIDPCRRep", fill="VariantID") {
  counts <- countsFlat %>%
            filter(SampleID %in% samples$SampleID) %>%
            filter(RefAllele == "True") %>%
            merge(samples %>% select("SampleID", group, "ControlForAmplicon","CellLine")) %>%
            dplyr:::rename(Frequency="%Reads") %>%
            mutate(Edited=ordered(ControlForAmplicon, levels=c(FALSE,TRUE), labels=c("Edited","Unedited"))) %>%
            as.data.frame()
  p <- ggplot(counts, aes_string(x=group, y="Frequency", fill=fill)) + geom_col(position=position_dodge())
  p <- p + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6)) +
        ylab("Ref Allele Frequency (%)")
  p <- p + theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
                 legend.title = element_text(size=9), #change legend title font size
                 legend.text = element_text(size=8)) #change legend text font size
  p <- p + facet_grid(cols=vars(Edited), scales = "free", space = "free")
  return(p)
}


pdf(file=paste0(opt$outbase, ".refAllele.barplots.pdf"), width=9, height=8)
samplesInput <- samplesheet %>% filter(Bin == "All" | ControlForAmplicon)
for (amplicon in unique(samplesheet$AmpliconID)) {
  currSamples <- samplesheet %>% filter(AmpliconID == amplicon)
  p <- getRefStackedBarplot(countsFlat, currSamples) + ggtitle(paste0("AmpliconID==",amplicon))
  print(p)
}
invisible(dev.off())



############################################
## Plot replicate concordance

flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
    )
}

getPCRReplicateCorrelations <- function(countsFlat, samplesheet, includeRef=FALSE) {
  ## Return correlations among all pairs of PCR replicates, considering editing rates of non-reference desired alleles
  results <- list()
  samplesheet <- samplesheet %>% mutate(Grouping=paste0(ExperimentIDReplicates))
  for (group in unique(samplesheet$Grouping)) {
    currSamples <- samplesheet %>% filter(Grouping == group)

    currCounts <- countsFlat %>%
      filter(SampleID %in% currSamples$SampleID) %>%
      filter(includeRef | RefAllele == "False") %>%
      select(SampleID,VariantID,`%Reads`) %>%
      spread(SampleID,`%Reads`,fill=0) %>%
      select(-VariantID) %>%
      as.matrix()

    correlations <- cor(currCounts)
    if (nrow(correlations) > 1) {
      cor.df <- flattenCorrMatrix(correlations)
      results[[group]] <- cor.df
    }
  }
  flat <- do.call(rbind, results)
  return(flat)
}


getPCRReplicateVariantCV <- function(countsFlat, samplesheet, includeRef=FALSE) {
  ## Return coefficient of variation for each variant among all pairs of PCR replicates — so that we can see the degree of variance as a function of allele frequency
  results <- list()
  samplesheet <- samplesheet %>% mutate(Grouping=paste0(ExperimentIDReplicates))
  for (group in unique(samplesheet$Grouping)) {
    currSamples <- samplesheet %>% filter(Grouping == group)

    if (nrow(currSamples) >= 2) {
      currCounts <- countsFlat %>%
        filter(SampleID %in% currSamples$SampleID) %>%
        filter(includeRef | RefAllele == "False") %>%
        select(SampleID,VariantID,`%Reads`) %>%
        spread(SampleID,`%Reads`,fill=0) %>%
        select(-VariantID) %>%
        as.matrix()

      cv <- t(apply(currCounts, 1, function(row) return(c(mean(row), sd(row), sd(row)/mean(row)*100))))
      results[[group]] <- cv
    }
  }
  flat <- data.frame(do.call(rbind, results))
  colnames(flat) <- c("mean", "sd", "CV")
  return(flat)
}


getReplicatePlot <- function(countsFlat, samplesheet) {
  cor.df <- getPCRReplicateCorrelations(countsFlat, samplesheet)
  vcv <- getPCRReplicateVariantCV(countsFlat, samplesheet)

  p1 <- ggplot(cor.df, aes(x=cor)) + geom_histogram(binwidth=0.01) +
      theme_classic() +
      xlab("Pearson correlation")

  p2 <- ggplot(vcv, aes(x=mean, y=CV)) + geom_point(alpha=0.5) +
      theme_classic() +
      xlab("Variant Frequency (%)") +
      ylab("Coefficient of Variation") +
      scale_x_log10()

  q <- plot_grid(p1, p2, ncol=2)
  return(q)
}

pdf(file=paste0(opt$outbase, ".replicateCorrelations.pdf"), width=6, height=3)
p <- getReplicatePlot(countsFlat, samplesheet)
print(p)
invisible(dev.off())



############################################
## Plot variant frequencies across bins

getBinnedBarplot <- function(countsFlat, samples, binList, group="ExperimentIDReplicates", normalizeBins=TRUE) {
  counts <- countsFlat %>%
            filter(SampleID %in% samples$SampleID) %>%
            merge(samples %>% select("SampleID",group,"Bin")) %>%
            filter(Bin %in% binList) %>%
            dplyr:::rename(Frequency="%Reads") %>%
            as.data.frame()

  countsGrouped <- counts %>%
            group_by_at(c(group, "VariantID", "Bin", "RefAllele")) %>%
            summarize(GroupedFrequency=mean(Frequency, na.rm=T)) %>%
            group_by_at(c(group,"VariantID","RefAllele")) %>%
            mutate(GroupedFrequencyAvgBin=mean(GroupedFrequency),
                   GroupedFreqRelativeToAverageBin=GroupedFrequency/mean(GroupedFrequency)) %>%
            as.data.frame()

  counts <- merge(counts, countsGrouped) %>% mutate(FreqRelativeToAverageBin=Frequency * GroupedFreqRelativeToAverageBin / GroupedFrequency)

  p <- countsGrouped %>% mutate(Variant=paste0(VariantID,"\n(",format(GroupedFrequencyAvgBin, digits=2),"%)")) %>%
       ggplot(aes(x=Variant, y=GroupedFreqRelativeToAverageBin, fill=Bin)) +
       geom_bar(stat="identity", position=position_dodge()) +
       scale_fill_grey(start=0.8, end=0.2) +
       #geom_point(data=counts, aes(x=VariantID, y=FreqRelativeToAverageBin, group=interaction(VariantID,Bin)), size=0.5, fill='red', position=position_dodge(width=0.25)) +
       geom_hline(yintercept=1, linetype="dashed", color="black") +
       ylim(0,1.5) + ylab("Frequency (normalized)") +
       theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6)) +
       theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
                 legend.title = element_text(size=9), #change legend title font size
                 legend.text = element_text(size=8)) #change legend text font size

  p2 <- counts %>% mutate(Variant=paste0(VariantID,"\n(",format(GroupedFrequencyAvgBin, digits=2),"%)")) %>%
       ggplot(aes_string(x="Variant", y="FreqRelativeToAverageBin", color="Bin")) +
       geom_boxplot(outlier.shape=NA, lwd=0.2) + geom_point(lwd=0.2, size=0.2, position=position_jitterdodge(jitter.width=0.05, seed=1)) +
       scale_fill_grey(start=0.8, end=0.2) + scale_color_grey(start=0.8, end=0.2) +
       geom_hline(yintercept=1, linetype="dashed", color="black") +
       ylim(0,1.5) + ylab("Frequency (normalized)") +
       theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6)) +
       theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
                 legend.title = element_text(size=9), #change legend title font size
                 legend.text = element_text(size=8)) #change legend text font size

  q <- plot_grid(p, p2, nrow=2)
  return(q)
}


pdf(file=paste0(opt$outbase, ".binBarplots.pdf"), width=7, height=8)
for (expt in unique(samplesheet$ExperimentIDReplicates)) {
  currSamples <- samplesheet %>% filter(ExperimentIDReplicates == expt)
  if (any(currSamples$Bin %in% binList)) {
    title <- ggdraw() + draw_label(paste0("ExperimentIDReplicates == ",expt), fontface='bold')
    p <- getBinnedBarplot(countsFlat, currSamples, binList, group="ExperimentIDReplicates")
    print(plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
  }
}
invisible(dev.off())


save.image(file=paste0(opt$outbase,".RData"))
