## HJ 3/2/2022
## Volcano Plot for every amplicon's variants
## Plot the power of the experiment (not partitioned by amplicon/target)


suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(stats4))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggExtra))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(pwr))


option.list <- list(
  make_option("--variantStatsFile", type="character", help="Variant stats file"),
  make_option("--volcanoPlots", type="character", help="Output PDF file with volcano plots"),
  make_option("--allEffects", type="character", help="File that contains all effects for every rep/experiment"),
  make_option("--powerPlots", type="character", help="Output PDF file with power plots"),
  make_option("--reps", type ="numeric", help="Combined number of BioReps and FFReps per variant")
)

opt <- parse_args(OptionParser(option_list=option.list))
if (is.null(opt$volcanoPlots))
  stop("PlotVariantStats: --volcanoPlots output should be specified.\n")

if (is.null(opt$powerPlots))
  stop("PlotVariantStats: --powerPlots output should be specified.\n")

if (is.null(opt$variantStatsFile) | !file.exists(opt$variantStatsFile))
  stop("PlotVariantStats: --variantStatsFile file not found.\n")

if (is.null(opt$allEffects) | !file.exists(opt$allEffects))
  stop("PlotVariantStats: --allEffects file not found.\n")


## Read data and add AmpliconInfo
VariantEffects <- read.delim(file = opt$variantStatsFile, sep = '\t')
allEffects <- read.delim(file = opt$allEffects, sep = '\t')
VariantEffects$AmpliconID <- allEffects$AmpliconID[match(VariantEffects[,1], allEffects$VariantID)]


## Local
# setwd('/Volumes/groups/engreitz/Projects/VariantEditing/FF/20230531-S019_P151_PPIF_enh_VFF/results/summary/')
# VariantEffects <- read.delim('AllelicEffectsStats.tsv', sep = '\t')
# allEffects <- read.delim(file = 'AllelicEffects.byExperimentRep.ExperimentIDReplicates.flat.tsv.gz', sep = '\t')
# VariantEffects$AmpliconID <- allEffects$AmpliconID[match(VariantEffects[,1], allEffects$VariantID)]


## Add 95% confidence interval and prepare for graphing
## Adding confidence intervals with specified number of reps
VariantEffects$Percent_Effect <- (VariantEffects$Mean_Effect - 1)*100
VariantEffects$Percent_SD <- VariantEffects$SD*100
VariantEffects$Percent_CI <- (1.96*(VariantEffects$SD/sqrt(4)))*100
VariantEffects$MinCells <- (VariantEffects$MinReps * VariantEffects$mean.cellCount)/4



## Add color scheme for differential effects of variants
VariantEffects$Expression_Effect <- "Insignificant"
VariantEffects$Expression_Effect[VariantEffects$Percent_Effect > 0.0 & VariantEffects$BH.p.value < 0.05] <- "Activating"
VariantEffects$Expression_Effect[VariantEffects$Percent_Effect < 0.0 & VariantEffects$BH.p.value < 0.05] <- "Suppressive"

## Volcano Plotting of Variant effects and significance
pdf(file=opt$volcanoPlots, width=9, height=7)
mycolors <- c("#00AFBB", "gray", "#E7B800")
names(mycolors) <- c("Activating", "Insignificant", "Suppressive")
for (amplicon in unique(VariantEffects$AmpliconID)) {
  ampliconEffects <- VariantEffects[VariantEffects$AmpliconID == amplicon, ]
  p1 <- ggplot(data=ampliconEffects, aes(x=Percent_Effect, y=-log10(BH.p.value), col=Expression_Effect)) +
    geom_point() +
    geom_errorbarh(aes(xmin=Percent_Effect-Percent_CI, xmax=Percent_Effect+Percent_CI)) +
    scale_color_manual(values = mycolors) +
    geom_hline(yintercept=-log10(0.05), col="red") +
    theme_minimal() +
    ylab("-Log10 Benjamini-Hochberg p-value") +
    xlab("% Effect on Gene Expression") +
    labs(col = "Effect") +
    ggtitle(label = paste(amplicon, "Variant Effects"))
  
  print(p1)
}

invisible(dev.off())



## Plotting necessary number of observations for every variant's detection
pdf(file=opt$powerPlots, width=8, height=8)
p2 <- ggplot(VariantEffects, aes(x=abs(Percent_Effect), y=MinCells)) + 
  geom_point(aes(col=Expression_Effect)) +
  geom_smooth(data = VariantEffects, method = 'loess', col = 'seashell4', fill = 'seashell2') +
  theme_minimal() +
  scale_y_continuous(trans='log10') +
  labs(col = "Effect") +
  scale_color_manual(values = mycolors) +
  xlab("Percent Effect on Gene Expression") +
  ylab("# of cells estimated to reach significance") +
  ggtitle(label = "Target number of cells per effect")

## Plotting estimated cell count vs observed cell count
high <- max(abs(VariantEffects$MinCells))
low <- min(abs(VariantEffects$MinCells))
p3 <- ggplot(VariantEffects, aes(x=mean.cellCount, y=MinCells, col=Expression_Effect)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, colour='black') +
  theme_minimal() +
  ylim(0,NA) +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  coord_fixed(xlim = c(low,high), ylim = c(low, high)) +
  labs(col = "Effect") +
  scale_color_manual(values = mycolors) +
  xlab("Number of cells in experiment") +
  ylab("Estimated target number of cells") +
  ggtitle(label = "Estimated target number of cells and the observed number of cells per variant")
  
## Plotting SD vs Effect size

VariantEffects$Proportion_SD <- VariantEffects$Percent_SD/VariantEffects$Percent_Effect
SigVars <- subset(VariantEffects, BH.p.value <0.05)
p4 <- ggplot(VariantEffects, aes(x=abs(Percent_Effect), y=abs(SD))) + 
  geom_point(aes(col=Expression_Effect)) +
  theme_minimal() +
  stat_cor() +
  #ylim(0, 10) + 
  #xlim(0, 20) +
  # geom_smooth(method = 'lm', se = FALSE, color = 'red', fullrange = TRUE) +
  labs(col = "Effect") +
  scale_color_manual(values = mycolors) +
  xlab("Percent Effect on Gene Expression") +
  ylab("SD")


p10 <- ggplot(SigVars, aes(x=abs(Percent_Effect), y=abs(SD))) + 
  geom_point(aes(col=Expression_Effect)) +
  theme_minimal() +
  stat_cor() +
  geom_smooth(method = 'lm', se = FALSE, color = 'red', fullrange = TRUE) +
  labs(col = "Effect") +
  scale_color_manual(values = mycolors) +
  xlab("Percent Effect on Gene Expression") +
  ylab("SD")

## SD vs cell number
p5 <- ggplot(SigVars, aes(x=freq, y=abs(SD))) + 
  geom_point(aes(col=Expression_Effect)) +
  theme_minimal() +
  scale_x_continuous(trans='log10') +
  stat_cor() +
  geom_smooth(method = 'lm', se = FALSE, color = 'red', fullrange = TRUE) +
  labs(col = "Effect") +
  scale_color_manual(values = mycolors) +
  xlab("Variant Frequency") +
  ylab("SD")

p11 <- ggplot(VariantEffects, aes(x=freq, y=abs(SD))) + 
  geom_point(aes(col=Expression_Effect)) +
  theme_minimal() +
  scale_x_continuous(trans='log10') +
  stat_cor() +
  geom_smooth(method = 'lm', se = FALSE, color = 'red', fullrange = TRUE) +
  labs(col = "Effect") +
  scale_color_manual(values = mycolors) +
  xlab("Variant Frequency") +
  ylab("SD")


## Power vs effect_size*number of cells
p6 <- ggplot(VariantEffects, aes(x=mean.cellCount*abs(Percent_Effect), y=Power)) + 
  geom_point(aes(col=Expression_Effect)) +
  ylim(0,1) + 
  theme_minimal() +
  geom_smooth(method = 'glm', col = 'seashell4', fill = 'seashell2', method.args = list(family = "binomial"), se = FALSE) +
  scale_color_manual(values = mycolors) +
  labs(col = "Effect") +
  xlab("% Effect * Number of cells") +
  ylab("Power")

## Plotting power vs estimated effect
p7 <- ggplot(VariantEffects, aes(x=abs(Percent_Effect), y=Power, col=Expression_Effect)) + 
  geom_point() +
  geom_errorbarh(aes(xmin=abs(Percent_Effect)-Percent_CI, xmax=abs(Percent_Effect)+Percent_CI)) +
  ylim(0,1) + 
  xlim(0,NA) +
  theme_minimal() +
  geom_smooth(method = 'glm', col = 'seashell4', fill = 'seashell2', method.args = list(family = "binomial"), se = FALSE) +
  scale_color_manual(values = mycolors) +
  labs(col = "Effect") +
  xlab("% Effect on Gene Expression") +
  ylab("Power")

## Plotting power vs number of cells observed
p8 <- ggplot(VariantEffects, aes(x=mean.cellCount, y=Power)) + 
  geom_point(aes(col=abs(Percent_Effect))) +
  scale_colour_gradient(low = 'white', high = 'red') +
  ylim(0,1) + 
  scale_x_continuous(trans='log10') +
  theme_minimal() +
  geom_smooth(method = 'glm', col = 'seashell4', fill = 'seashell2', method.args = list(family = "binomial"), se = FALSE) +
  labs(col = "Estimated % Effect on Gene Expression") +
  xlab("Number of cells per FFRep") +
  ylab("Power")


## Plotting Power vs cell # w/o gradient color
p9 <- ggplot(VariantEffects, aes(x=mean.cellCount, y=Power)) + 
  geom_point(aes(col=Expression_Effect)) +
  ylim(0,1) + 
  scale_x_continuous(trans='log10') +
  theme_minimal() +
  scale_color_manual(values = mycolors) +
  geom_smooth(method = 'glm', col = 'seashell4', fill = 'seashell2', method.args = list(family = "binomial"), se = FALSE) +
  labs(col = "Effect") +
  xlab("Number of cells per FFRep") +
  ylab("Power")


## Histogram of Variant Counts
# p11 <- ggplot(VariantEffects, aes(x = mean.cellCount)) + 
#   geom_histogram() + 
#   scale_x_log10() + 
#   theme_minimal() +
#   xlab("Cells with Variant per FFRep")


print(p4)
print(p5)
print(p2)
print(p3)
print(p6)
print(p7)
print(p8)
print(p9)
print(p10)
print(p11)
invisible(dev.off())






