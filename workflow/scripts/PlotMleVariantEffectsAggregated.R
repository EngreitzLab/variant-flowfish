# Jesse Engreitz, Katherine Guo
# 3/10/22
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
    opt$allelicEffectFile = "results/summary/AllelicEffects.byPCRRep.ExperimentIDPCRRep.flat.tsv.gz"
    opt$samplesheet = "SampleList.snakemake.tsv"
    opt$outfile = "results/summary/AllelicEffects.byPCRRep.ExperimentIDPCRRep.pdf"
    opt$effectColumn = "effect_size"
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
suppressPackageStartupMessages(library(Hmisc))

mle <- read.delim(opt$allelicEffectFile, check.names=F, stringsAsFactors=F)
stopifnot(opt$effectColumn %in% colnames(mle))
mle <- subset(mle, sum1 >= 1000) ## Only plot the variants with at least 1000 cells

samplesheet <- read.delim(opt$samplesheet, check.names=F, stringsAsFactors=F)

pdf(file=opt$outfile, width=12, height=6)

for (amplicon in unique(samplesheet$AmpliconID)) {
    p1 <- ggplot(mle %>% filter(AmpliconID==amplicon), aes_string("VariantID", opt$effectColumn)) +
        geom_point(size=0.05) + 
        geom_jitter(width = 0.05) +
        stat_summary(aes(VariantID, effect_size), fun.data = mean_cl_normal, geom="errorbar") + 
        stat_summary(fun.y=mean, geom="point", shape=18,
                        size=2, color="red") + 
        ylab("Effect size on gene expression\n(MLE, vs reference allele)") +
        ggtitle(paste(amplicon, "Variant Effects")) +
        geom_hline(yintercept=1, linetype="dashed", color="red") + 
        scale_y_continuous(limits=c(0, 1.5)) + 
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6), plot.title = 
        element_text(hjust = 0.5)) +
        theme(strip.background = element_blank(), strip.text.x = element_blank()) # don't show facet grid labels
    if (length(unique((mle %>% filter(AmpliconID==amplicon))$VariantID)) > 10) {
        p1 <- p1 + facet_grid(cols = vars(Location), scales = "free_x", switch = "x", space = "free_x")
    }
    print(p1)
  
    p2 <- mle %>% filter(AmpliconID==amplicon) %>% mutate(PctEffect=(get(opt$effectColumn)-1)*100) %>%
        ggplot(aes_string(x="VariantID", y="PctEffect")) +
        geom_point(size=0.05) + 
        geom_jitter(width = 0.05) +
        stat_summary(aes(VariantID, PctEffect), fun.data = mean_cl_normal, geom="errorbar") + 
        stat_summary(fun.y=mean, geom="point", shape=18,
                        size=2, color="red") + 
        ylab("Effect size on gene expression\n(MLE, % change vs reference allele)") +
        ggtitle(paste(amplicon, "Variant Effect Percent Change")) +
        geom_hline(yintercept=0, linetype="dashed", color="red") + 
        # facet_grid(cols = vars(Location), scales = "free_x", switch = "x", space = "free_x") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6), plot.title = 
        element_text(hjust = 0.5)) +
        theme(strip.background = element_blank(), strip.text.x = element_blank()) # don't show facet grid labels
    if (length(unique((mle %>% filter(AmpliconID==amplicon))$VariantID)) > 10) {
        p2 <- p2 + facet_grid(cols = vars(Location), scales = "free_x", switch = "x", space = "free_x")
    }
    print(p2)
}
invisible(dev.off())

save.image(file=paste0(opt$outfile,".RData"))
