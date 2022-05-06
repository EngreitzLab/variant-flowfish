## Hank Jones
## 3/2/2022
## Code for performing one sample two-tailed T-test and for estimating sample power
## Control for multiple testing using Benjamini-Hochberg corrected p-values

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(pwr))
suppressPackageStartupMessages(library(ggExtra))
suppressPackageStartupMessages(library(ggpubr))



option.list <- list(
  make_option("--allelicEffectFile", type="character", help="Allelic effect flat file"),
  make_option("--outfile", type="character", help="Output stats table"),
  make_option("--effectColumn", type="character", default="effect_size"),
  make_option("--repPlots", type="character", help="Output file name for replicate plots")
)
opt <- parse_args(OptionParser(option_list=option.list))

if (is.null(opt$outfile))
  stop("t.test: --outfile should be specified.\n")

if (is.null(opt$allelicEffectFile) | !file.exists(opt$allelicEffectFile))
  stop("t.test: --allelicEffectFile file not found.\n")


## Read data 
mle <- read.delim(opt$allelicEffectFile, check.names=F, stringsAsFactors=F)
stopifnot(opt$effectColumn %in% colnames(mle))
mle$Percent_Effect <- (mle$effect_size - 1)*100


## Drop references and variants without observations
mle <- subset(mle, effect_size != 1 & logMean != 0 & sum1 > 500)

## Capture the highest observed value to set limits of graphs
high <- max(abs(mle$Percent_Effect))

## Plot Replicates
pdf(file=opt$repPlots, width=9, height=9)

corData <- unique.data.frame(mle['VariantID'])

for (rep in unique(mle$BioRep)) {
  counter <- paste0("BioRep", rep)
  rep <- mle[mle$BioRep == rep, c('VariantID', 'Percent_Effect')]
  rep <- aggregate(Percent_Effect ~ ., mean, data = rep)
  corData <- merge(corData, rep)
  corData <- corData %>% rename(!!counter := Percent_Effect)
  
}

p1 <- ggplot(corData, aes(x=BioRep1, y=BioRep2)) +
  geom_point() +
  scale_x_continuous(expand=c(0,0), limits=c(-high-(high/10), high+(high/10))) +
  scale_y_continuous(expand=c(0,0), limits=c(-high-(high/10), high+(high/10))) +
  coord_cartesian(xlim=c(-high-(high/10), high+(high/10)), ylim=c(-high-(high/10), high+(high/10))) +
  geom_smooth(method = 'lm', se = FALSE, color = 'red', fullrange = TRUE) +
  geom_abline(intercept = 0, slope = 1, colour='black') +
  geom_hline(yintercept = 0, colour='seashell4') +
  geom_vline(xintercept = 0, colour='seashell4') +
  theme_minimal() +
  stat_cor() +
  xlab("BioRep 1 % Effect") +
  ylab("BioRep 2 % Effect") + 
  ggtitle(label = "BioReplicate Correlation")

p2 <- ggMarginal(p1, type = 'density', color = 'seashell3', fill = adjustcolor( "seashell3", alpha.f = 0.2))
print(p2)


for (rep in unique(mle$BioRep)) {
  FF1 <- mle[mle$BioRep == rep & mle$FFRep == 1, c('VariantID', 'Percent_Effect')]
  FF2 <- mle[mle$BioRep == rep & mle$FFRep == 2, c('VariantID', 'Percent_Effect')]
  corData <- merge(FF1, FF2, by = 'VariantID', all = TRUE)
  corData[is.na(corData)] <- 1
  p3 <- ggplot(corData, aes(x=Percent_Effect.x, y=Percent_Effect.y)) +
    geom_point() +
    stat_cor() +
    scale_x_continuous(expand=c(0,0), limits=c(-high-(high/10), high+(high/10))) +
    scale_y_continuous(expand=c(0,0), limits=c(-high-(high/10), high+(high/10))) +
    coord_cartesian(xlim=c(-high-(high/10), high+(high/10)), ylim=c(-high-(high/10), high+(high/10))) +
    geom_smooth(method = 'lm', se = FALSE, color = 'red', fullrange = TRUE) +
    geom_abline(intercept = 0, slope = 1, colour='gray1') +
    geom_hline(yintercept = 0, colour='seashell4') +
    geom_vline(xintercept = 0, colour='seashell4') +
    theme_minimal() +
    xlab("FlowFISH Rep1 % Effect") +
    ylab("FlowFISH Rep2 % Effect") + 
    ggtitle(label = paste0("FFReplicate Correlation: BioRep", rep))

  print(p3)
}

invisible(dev.off())

## Compute T-test for every variant and mean effect size
## Build a new dataframe as an output
## Compute number of observations needed at power of 0.9
## Use Bonferroni corrected p-value

variantStats <- data.frame(row.names = unique(mle$VariantID))
test.p.val <- 0.05/nrow(variantStats)

for (variant in unique(mle$VariantID)) {
  t.data <- mle[(mle$VariantID == variant & mle$sum1 > 500),]
  t.data <- t.data[, c('AmpliconID', 'sum1', 'effect_size')]
  p.value <- 1
  MinReps <- NA
  Power <- NA
  SD <- sd(t.data$effect_size)
  Mean_Effect <- mean(t.data$effect_size)
  mean.cellCount <- mean(t.data$sum1)
  sd.cellCount <- sd(t.data$sum1)
  total.cellCount <- as.integer(sum(t.data$sum1))
  variantDelta <- (abs(1-mean(t.data$effect_size)))/SD
  if (length(t.data$effect_size) > 1) {
    if (mean(t.data$effect_size) != 1){
      tmp <- t.test(t.data$effect_size, mu = 1)
      p.value <- tmp$p.value
      pwr.tmp <- pwr.t.test(d = variantDelta, sig.level = test.p.val, power = 0.9, type = "one.sample")
      MinReps <- pwr.tmp$n
      pwr.tmp2 <- pwr.t.test(d = variantDelta, sig.level = test.p.val, n = (total.cellCount/mean.cellCount), type = 'one.sample')
      Power <- pwr.tmp2$power
    }
  }
  variantStats[variant, c('Mean_Effect', 'p.value', 'SD', 'MinReps', 'Power', 'mean.cellCount', 'sd.cellCount', 'total.cellCount')] <- c(Mean_Effect, p.value, SD, MinReps, Power, mean.cellCount, sd.cellCount, total.cellCount)
}


## Correct p values with BH
variantStats["BH.p.value"] <- p.adjust(variantStats$p.value, "fdr", nrow(variantStats))


## Write data
write.table(variantStats, file = opt$outfile, sep='\t', col.names = NA)
save.image(file=paste0(opt$outfile,".RData"))
