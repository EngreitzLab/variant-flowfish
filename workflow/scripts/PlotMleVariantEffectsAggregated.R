# Jesse Engreitz
# 12/15/21
# Rscript to plot MLE effect size estimates, combined across samples


suppressPackageStartupMessages(library("optparse"))

option.list <- list(
  make_option("--allelicEffectFile", type="character", help="Allelic effect flat file"),
  make_option("--samplesheet", type="character", help="Snakemake sample sheet including SampleID info"),
  make_option("--outfile", type="character", help="Output PDF file with plots"),
  make_option("--effectColumn", type="character", default="effect_size")
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


if (is.null(opt$outfile))
  stop("PlotMleVariantEffectsAggregated: --outfile should be specified.\n")

if (is.null(opt$allelicEffectFile) | !file.exists(opt$allelicEffectFile))
  stop("PlotMleVariantEffectsAggregated: --allelicEffectFile file not found.\n")

if (is.null(opt$samplesheet) | !file.exists(opt$samplesheet))
    stop("PlotMleVariantEffectsAggregated: --samplesheet file not found.\n")

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))


mle <- read.delim(opt$allelicEffectFile, check.names=F, stringsAsFactors=F)
stopifnot(opt$effectColumn %in% colnames(mle))

samplesheet <- read.delim(opt$samplesheet, check.names=F, stringsAsFactors=F)


p1 <- ggplot(mle, aes_string(x="VariantID", y=opt$effectColumn)) +
     geom_col() + 
     ylab("Effect size on gene expression\n(MLE, vs reference allele)") +
     geom_hline(yintercept=1, linetype="dashed", color="black") + 
     theme_classic() +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6)) 

p2 <- mle %>% mutate(PctEffect=(get(opt$effectColumn)-1)*100) %>%
     ggplot(aes_string(x="VariantID", y="PctEffect")) +
     geom_col() + 
     ylab("Effect size on gene expression\n(MLE, % change vs reference allele)") +
     geom_hline(yintercept=0, linetype="dashed", color="black") + 
     theme_classic() +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6)) 

pdf(file=opt$outfile, width=3, height=4)
print(p1)
print(p2)
dev.off()






############################################
## Plot variant frequencies across bins

getBinnedBarplot <- function(mle, samples, idcol) {

  counts <- mle %>%
            merge(samples %>% select_at(c(idcol,"AmpliconID")) %>% unique()) %>%
            group_by(VariantID) %>%
            mutate(GroupedFrequencyAvg=mean(freq),
                   Variant=paste0(VariantID,"\n(",format(GroupedFrequencyAvg, digits=2, scientific=FALSE),"%)"),
                   Replicate=factor(1:n()),
                   PctEffect=(get(opt$effectColumn)-1)*100) %>%
            as.data.frame()

  p <- counts %>% 
       ggplot(aes_string(x="Variant", y="PctEffect", fill="Replicate")) +
       geom_bar(stat="identity", position=position_dodge()) +
       scale_fill_grey(start=0.8, end=0.2) +
       geom_hline(yintercept=0, linetype="dashed", color="black") +
       ylab("Effect Size (% vs reference allele)") +
       theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6)) +
       theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
                 legend.title = element_text(size=9), #change legend title font size
                 legend.text = element_text(size=8)) #change legend text font size

  p2 <- counts %>%
       ggplot(aes_string(x="Variant", y="PctEffect")) +
       geom_boxplot(outlier.shape=NA, lwd=0.2) + geom_point(lwd=0.2, size=0.2) + 
       scale_fill_grey(start=0.8, end=0.2) + scale_color_grey(start=0.8, end=0.2) +
       geom_hline(yintercept=0, linetype="dashed", color="black") +
       ylab("Effect Size (% vs reference allele)") +
       theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6)) +
       theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
                 legend.title = element_text(size=9), #change legend title font size
                 legend.text = element_text(size=8)) #change legend text font size

  q <- plot_grid(p, p2, nrow=2)
  return(q)
}


idcol <- intersect(colnames(mle), colnames(samplesheet)) %>% setdiff("AmpliconID")

pdf(file=paste0(opt$outfile), width=7, height=8)
for (expt in unique(samplesheet$ExperimentID)) {
  currSamples <- samplesheet %>% filter(ExperimentID == expt)

  if (any(currSamples[,idcol] %in% mle[,idcol])) {
    title <- ggdraw() + draw_label(paste0("ExperimentID == ",expt), fontface='bold')
    p <- getBinnedBarplot(mle, currSamples, idcol)
    print(plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
  }
}
invisible(dev.off())


save.image(file=paste0(opt$outfile,".RData"))
