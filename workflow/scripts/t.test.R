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
  make_option("--varPlots", type="character", help = "Name of variance modeling/power calculation plots"),
  make_option("--volcanoPlots", type="character", help="Output PDF file with volcano plots")
  # make_option("--repPlots", type="character", help="Output file name for replicate plots")
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

## Local run
# setwd("/Users/meat/Downloads/")
# mle <- read.delim("AllelicEffects.byExperimentRep.ExperimentIDReplicates.flat.tsv.gz", check.names = FALSE, stringsAsFactors =F)
# mle$Percent_Effect <- (mle$effect_size - 1)*100

## Drop references and variants without observations below freq 0.0001
mle <- subset(mle, effect_size != 1 & logMean != 0 & freq > 0.0001)

## Capture the highest observed value to set limits of graphs
high <- max(abs(mle$Percent_Effect))

## Plot Replicates
# pdf(file=opt$repPlots, width=9, height=9)
# 
# corData <- unique.data.frame(mle['VariantID'])
# 
# for (rep in unique(mle$BioRep)) {
#   counter <- paste0("BioRep", rep)
#   rep <- mle[mle$BioRep == rep, c('VariantID', 'Percent_Effect')]
#   rep <- aggregate(Percent_Effect ~ ., mean, data = rep)
#   corData <- merge(corData, rep)
#   corData <- corData %>% rename(!!counter := Percent_Effect)
#   
# }
# 
# p1 <- ggplot(corData, aes(x=BioRep1, y=BioRep2)) +
#   geom_point() +
#   scale_x_continuous(expand=c(0,0), limits=c(-high-(high/10), high+(high/10))) +
#   scale_y_continuous(expand=c(0,0), limits=c(-high-(high/10), high+(high/10))) +
#   coord_cartesian(xlim=c(-high-(high/10), high+(high/10)), ylim=c(-high-(high/10), high+(high/10))) +
#   geom_smooth(method = 'lm', se = FALSE, color = 'red', fullrange = TRUE) +
#   geom_abline(intercept = 0, slope = 1, colour='black') +
#   geom_hline(yintercept = 0, colour='seashell4') +
#   geom_vline(xintercept = 0, colour='seashell4') +
#   theme_minimal() +
#   stat_cor() +
#   xlab("BioRep 1 % Effect") +
#   ylab("BioRep 2 % Effect") + 
#   ggtitle(label = "BioReplicate Correlation")
# 
# p2 <- ggMarginal(p1, type = 'density', color = 'seashell3', fill = adjustcolor( "seashell3", alpha.f = 0.2))
# print(p2)
# 
# 
# for (rep in unique(mle$BioRep)) {
#   FF1 <- mle[mle$BioRep == rep & mle$FFRep == 1, c('VariantID', 'Percent_Effect')]
#   FF2 <- mle[mle$BioRep == rep & mle$FFRep == 2, c('VariantID', 'Percent_Effect')]
#   corData <- merge(FF1, FF2, by = 'VariantID', all = TRUE)
#   corData[is.na(corData)] <- 1
#   p3 <- ggplot(corData, aes(x=Percent_Effect.x, y=Percent_Effect.y)) +
#     geom_point() +
#     stat_cor() +
#     scale_x_continuous(expand=c(0,0), limits=c(-high-(high/10), high+(high/10))) +
#     scale_y_continuous(expand=c(0,0), limits=c(-high-(high/10), high+(high/10))) +
#     coord_cartesian(xlim=c(-high-(high/10), high+(high/10)), ylim=c(-high-(high/10), high+(high/10))) +
#     geom_smooth(method = 'lm', se = FALSE, color = 'red', fullrange = TRUE) +
#     geom_abline(intercept = 0, slope = 1, colour='gray1') +
#     geom_hline(yintercept = 0, colour='seashell4') +
#     geom_vline(xintercept = 0, colour='seashell4') +
#     theme_minimal() +
#     xlab("FlowFISH Rep1 % Effect") +
#     ylab("FlowFISH Rep2 % Effect") + 
#     ggtitle(label = paste0("FFReplicate Correlation: BioRep", rep))
# 
#   print(p3)
# }
# 
# invisible(dev.off())

## Compute T-test for every variant and mean effect size
## Build a new dataframe as an output
## Compute number of observations needed at power of 0.9
## Use Bonferroni corrected p-value

variantStats <- data.frame(row.names = unique(mle$VariantID))
test.p.val <- 0.05/nrow(variantStats)

for (variant in unique(mle$VariantID)) {
  ## Select variants with at least 500 alleles per rep
  ## Choose the AmpliconID, total number of cells in rep, and the estimated effect size
  t.data <- mle[(mle$VariantID == variant & mle$freq > 0.0001),]
  t.data <- t.data[, c('AmpliconID', 'sum1', 'effect_size', 'freq')]
  ## Initialize values for power code
  ## Calculate statistics
  reps <- length(t.data$AmpliconID)
  p.value <- 1
  MinReps <- NA
  Power <- NA
  SD <- sd(t.data$effect_size)
  Variance <- var(t.data$effect_size)
  freq <- mean(t.data$freq)
  Mean_Effect <- mean(t.data$effect_size)
  mean.cellCount <- mean(t.data$sum1)
  sd.cellCount <- sd(t.data$sum1)
  total.cellCount <- as.integer(sum(t.data$sum1))
  ## Calculate Cohen's D
  variantDelta <- (abs(1-mean(t.data$effect_size)))/SD
  ## Select only variants with observations in more than 1 rep
  if (length(t.data$effect_size) > 1) {
    ## Get rid of reference alleles
    if (mean(t.data$effect_size) != 1){
      ## Run one sample t-test for variant.
      ## Null is x = 1
      if (length(unique(t.data$effect_size)) != 1){ # error when all of the effects are the same
        tmp <- t.test(t.data$effect_size, mu = 1)
        p.value <- tmp$p.value
        ## Estimate power
        ## Use the estimated effect size, sig.level is the bonferroni corrected p-val, power at 0.9
        ## Power at 0.9 means that 90% of the time the negatives are true negatives
        ## Compute sample size needed to pass significance (ie MinReps)
        pwr.tmp <- pwr.t.test(d = variantDelta, sig.level = test.p.val, power = 0.9, type = "one.sample")
        MinReps <- pwr.tmp$n
        ## Estimate power using effect sizes and standard deviations
        pwr.tmp2 <- pwr.t.test(d = variantDelta, sig.level = test.p.val, n = (total.cellCount/mean.cellCount), type = 'one.sample')
        Power <- pwr.tmp2$power
      }
    }
  }
  variantStats[variant, c('Mean_Effect', 'p.value', 'SD', 'Variance', 'MinReps', 'Power', 'mean.cellCount', 'sd.cellCount', 'freq', 'total.cellCount', 'reps')] <- c(Mean_Effect, p.value, SD, Variance, MinReps, Power, mean.cellCount, sd.cellCount, freq, total.cellCount, reps)
}


## Correct p values with BH
variantStats["BH.p.value"] <- p.adjust(variantStats$p.value, "fdr", nrow(variantStats))


## Write data
write.table(variantStats, file = opt$outfile, sep='\t', col.names = NA)
save.image(file=paste0(opt$outfile,".RData"))

## Local
## write.table(variantStats, file = "AllelicEffectsStats.tsv", sep='\t', col.names = NA)

################################################################################################

## Modeling variance as a function of variant frequency


library(gridExtra)
variantStats$Percent_Effect <- (variantStats$Mean_Effect - 1)*100
variantStats <- na.omit(variantStats)
variantStats$abs_Effect <- abs(variantStats$Percent_Effect)
variantStats$Percent_Variance <- variantStats$Variance * 100

model <- lm(log(Percent_Variance) ~ log(freq), data = variantStats)

# Generate new data for predictions
newdata <- data.frame(
  freq = exp(seq(min(log(variantStats$freq)), max(log(variantStats$freq)), length.out = 100))
)

# Add predicted log(Percent_Variance) to the newdata
newdata$predicted_log_variance <- predict(model, newdata = newdata)

# Convert predicted log(Percent_Variance) back to the original scale
newdata$predicted_Percent_Variance <- exp(newdata$predicted_log_variance)

# Plotting
mytheme <- theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))

p10 <- ggplot(variantStats, aes(x = freq, y = Percent_Variance)) + 
  geom_point() + 
  geom_line(data = newdata, aes(x = freq, y = predicted_Percent_Variance), color = "red") +
  labs(title = "Variance vs Frequency") +
  xlab("Variant Frequency") +
  ylab("Variance") +
  xlim(0, 0.02) + 
  ylim(0, 8)
p10 <- p10 + mytheme
p10

## Log Space plotting
newdata$log_freq = log(newdata$freq)

p11 <- ggplot(variantStats, aes(x = log(freq), y = log(Percent_Variance))) + 
  geom_point() + 
  geom_line(data = newdata, aes(x = log_freq, y = predicted_log_variance), color = "red") +
  labs(title = "Log(Variance) vs Log(Frequency)") +
  xlab("Log(Variant Frequency)") +
  ylab("Log(Variance)")
p11 <- p11 + mytheme
p11


#######################################################################################

## Computing power at various effect sizes as a function of variant frequencies

compute_power <- function(effect_size, n_replicates, freq) {
  predicted_variance <- predict(model, newdata = data.frame(freq = freq))
  sd <- sqrt(exp(predicted_variance)) # transforming variance back to original scale
  power <- power.t.test(n = n_replicates, delta = effect_size, sd = sd, 
                        sig.level = 0.05, type = "two.sample", alternative = "two.sided")$power
  return(power)
}

# Compute power for effect sizes: 8 replicates = number of FFReps
effect_sizes <- seq(0.1, 0.5, by = 0.05)
n_replicates <- 8
freq <- frequencies <- seq(0.001, 0.05, by = 0.001)

# Initialize a data frame to store the results
power_df <- data.frame()

# Loop over frequencies and effect sizes
for (f in frequencies) {
  for (e in effect_sizes) {
    power <- compute_power(effect_size = e, n_replicates = n_replicates, freq = f)
    power_df <- rbind(power_df, data.frame(Frequency = f, EffectSize = e, Power = power))
  }
}

# Plotting
p12 <- ggplot(power_df, aes(x = Frequency, y = Power)) +
  geom_line() +
  facet_wrap(~EffectSize, scales = "free") +
  xlab("Frequency") +
  ylab("Power") +
  theme_minimal() +
  ggtitle("Power as a function of Frequency and Effect Size")

p12 <- p12 + mytheme
p12




pdf(file=opt$varPlots, width=9, height=9)
p10
p11
p12
invisible(dev.off())


## Volcano Plots
variantStats$AmpliconID <- mle$AmpliconID[match(rownames(variantStats), mle$VariantID)]

## Add 95% confidence interval and prepare for graphing
## Adding confidence intervals with specified number of reps
variantStats$Percent_Effect <- (variantStats$Mean_Effect - 1)*100
variantStats$Percent_SD <- variantStats$SD*100
variantStats$Percent_CI <- (1.96*(variantStats$SD/sqrt(variantStats$reps)))*100
variantStats$MinCells <- (variantStats$MinReps * variantStats$mean.cellCount)/variantStats$reps



## Add color scheme for differential effects of variants
variantStats$Expression_Effect <- "Insignificant"
variantStats$Expression_Effect[variantStats$Percent_Effect > 0.0 & variantStats$BH.p.value < 0.05] <- "Activating"
variantStats$Expression_Effect[variantStats$Percent_Effect < 0.0 & variantStats$BH.p.value < 0.05] <- "Suppressive"

## Volcano Plotting of Variant effects and significance
pdf(file=opt$volcanoPlots, width=9, height=7)
mycolors <- c("#00AFBB", "gray", "#E7B800")
names(mycolors) <- c("Activating", "Insignificant", "Suppressive")
for (amplicon in unique(variantStats$AmpliconID)) {
  ampliconEffects <- variantStats[variantStats$AmpliconID == amplicon, ]
  p13 <- ggplot(data=ampliconEffects, aes(x=Percent_Effect, y=-log10(BH.p.value), col=Expression_Effect)) +
    geom_point() +
    geom_errorbarh(aes(xmin=Percent_Effect-Percent_CI, xmax=Percent_Effect+Percent_CI)) +
    scale_color_manual(values = mycolors) +
    geom_hline(yintercept=-log10(0.05), col="red") +
    theme_minimal() +
    ylab("-Log10 Benjamini-Hochberg p-value") +
    xlab("% Effect on Gene Expression") +
    labs(col = "Effect") +
    ggtitle(label = paste(amplicon, "Variant Effects"))
}
p13 <- p13 + mytheme
p13
invisible(dev.off())










