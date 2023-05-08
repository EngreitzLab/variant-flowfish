# Jesse Engreitz
# 5/30/21
# Rscript to plot MLE effect size estimates


suppressPackageStartupMessages(library("optparse"))

option.list <- list(
  make_option("--mleEffects", type="character", help="MLE effect sizes"),
  make_option("--effectColumn", type="character", default="effect_size"),
  make_option("--outfile", type="character", help="Output plot filename")
  )
opt <- parse_args(OptionParser(option_list=option.list))
dput(opt)


suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))


mle <- read.delim(opt$mleEffects, check.names=F, stringsAsFactors=F)
mle <- subset(mle, VariantID != "")  ## Only plot results for desired variants
mle <- subset(mle, sum1 >= 1000) ## Only plot the variants with at least 1000 cells

stopifnot(opt$effectColumn %in% colnames(mle))

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
