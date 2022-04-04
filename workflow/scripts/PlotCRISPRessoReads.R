options(scipen=999)
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))

option.list <- list(
  make_option("--CRISPRessoAggFolder", type="character", help="CRISPResso Aggregate Directory"),
  make_option("--genotyping_only", type="logical", help="If analysis is genotyping only"),
  make_option("--outfile", type="character", help="Output plot filename")
  )
opt <- parse_args(OptionParser(option_list=option.list))

crispresso_agg <- paste0(opt$CRISPRessoAggFolder, "/CRISPRessoAggregate_quantification_of_editing_frequency_by_amplicon.txt")

agg <- read.delim(crispresso_agg, check.names=F, stringsAsFactors=F)

if (opt$genotyping_only) {
  agg <- extract(agg, Folder, into = c("Sample"), "^.\\/CRISPResso_on_(.*)")
} else {
  agg <- extract(agg, Folder, into = c("Sample", "Bin"), "^.\\/CRISPResso_on_(.*)-(Bin.*$)")
}


pdf(file=opt$outfile, width=12, height=6)
p <- agg %>% gather("Reads", "num_reads", c('Reads_in_input', 'Reads_aligned')) %>% 
    ggplot(aes(x = Sample,y = num_reads)) + 
    geom_bar(aes(fill = Reads),stat = "identity",position = "dodge") + 
    ggtitle("CRISPResso Aggregate Read Counts") +
    theme_classic() +
    xlab("Sample ID") +
    ylab("# Reads") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6), plot.title = 
        element_text(hjust = 0.5)) 
print(p)
invisible(dev.off())
