#!/usr/bin/env Rscript
##################################################################
## EstimateEffectsFromBins.R
## Jesse Engreitz and Ben Doughty
##
## Based on previous script from Tejal Patwardhan, with a critical fix to the log likelihood function
## July 7, 2017
## 
## Input: count data, sort params, design file
## Black box: computes weighted average, mle mean, standard deviation for each guides
## Output: .bed and .bedgraph for each guide that passes filters
##

rm(list=ls())
suppressPackageStartupMessages(library("optparse"))
option.list <- list(
  make_option(c("-dlo", "--designDocLocation"), type="character", help="Design document location, ex: ~/Documents/Harvard/LanderResearch/flo1/raw/160705.Ess.Design.txt"),
  make_option(c("-clo", "--countsLocation"), type="character", help="Counts document location, ex: ~/Documents/Harvard/LanderResearch/flo1/raw/4.tsv"),
  make_option(c("-sp", "--sortParamsloc"), type="character", help="MUST HAVE A MEAN, BOUNDS, AND BARCODE COLUMN!!, ex: ~/Documents/Harvard/LanderResearch/flo1/raw/sortParams/gata1ff3.txt"),
  make_option(c("-om", "--outputmle"), type="character", help="File to write the MLE mean into for each guide"),
  make_option(c("-l", "--log"), type="character", help="Log results of MLE")
)
opt <- parse_args(OptionParser(option_list=option.list))

## Set up variables
designDocLocation <- opt$designDocLocation
countsLocation <- opt$countsLocation
sortParamsloc <- opt$sortParamsloc
outputmle <- opt$outputmle
log <- opt$log

MAXMEAN <- 10000 # in FACS space, max value a guide can take
MINMEAN <- 1  # in FACS space, max value a guide can take

## load packages
# suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stats4))
suppressPackageStartupMessages(library(methods))
set.seed(100)

## FUNCTIONS

## NLL function
# takes as input a mean and standard deviation in log space, as well as a set of bin boundaries and a set of observations
# per bin (which is of length nbins+1, to account for missed cells)
# returns the negative of the log likelihood of observing those counts with the given parameters
ll <- function(mu, std, observations, bins) {
  write("Invoked NLL",file=log,append=TRUE)
  # write(dim(bins),file=log,append=TRUE)

  #writeLines("Invoked NLL", log)
  if (std < 0) {
    return(10^10) ## make the sum really high if it picks a negative standard deviation
  }
  
  pe <- vector()
  
  for (i in c(1:dim(bins)[1])) {
    write(i,file=log,append=TRUE)
    pe <- c(pe, pnorm(bins[i,2], mean=mu, sd=std, log.p=FALSE) - pnorm(bins[i,1], mean=mu, sd=std, log.p=FALSE))
  }  

  ## Probabilities of each bin, as determined by a CDF function for a normal distribution
  #pe <- c(pnorm(bins[1,2], mean=mu, sd=std, log.p=FALSE),
  #        pnorm(bins[2,2], mean=mu, sd=std, log.p=FALSE) - pnorm(bins[2,1], mean=mu, sd=std, log.p=FALSE),
  #        pnorm(bins[3,2], mean=mu, sd=std, log.p=FALSE) - pnorm(bins[3,1], mean=mu, sd=std, log.p=FALSE),
  #        pnorm(bins[4,2], mean=mu, sd=std, log.p=FALSE) - pnorm(bins[4,1], mean=mu, sd=std, log.p=FALSE),
  #        pnorm(bins[5,2], mean=mu, sd=std, log.p=FALSE) - pnorm(bins[5,1], mean=mu, sd=std, log.p=FALSE),
  #        pnorm(bins[6,1], mean=mu, sd=std, log.p=FALSE, lower.tail=FALSE))
  pe <- c(pe, 1-sum(pe)) # add a "bin" for the remaining cells

  pe[pe == 0] <- 10^-10 # remove any 0s from the probabilities (shouldn't happen but might)

  sum <- 0  # -log likelihood function
  for (i in c(1:length(observations))) {
    sum <- sum - (log(pe[i]) * observations[i])
  }

  return(sum)
}

# Wrap the maximum likelihood estimator
getNormalMLE <- function(mu.i, sd.i, bin.counts, total.count, bins, minmean=MINMEAN, maxmean=MAXMEAN, minvar=0, maxvar=1) {
  ## mu.i = initial guess of the mean
  ## bin.counts = vector of counts per bin
  ## total.count = estimated number of total counts (including cells not included in any bin)
  ## bins = matrix of bin boundaries (nrows = nbins, ncols = 2)
  ## returns mean and standard deviation in log-space
  bin.counts <- as.numeric(as.matrix(bin.counts))

  if (TRUE){#total.count < sum(bin.counts)) {
    o <- c(bin.counts, 0)
    #print(o)

    est <- mle(minuslog=ll, 
             start=list(mu=mu.i,std=sd.i), 
             fixed=list(observations=o, bins=bins), 
             method="L-BFGS-B",
             lower=c(log10(minmean), minvar),
             # lower=c(0, 0), since we are in log space, I think there are no constraints on the mean 
             upper=c(log10(maxmean), maxvar))#, 
             #control=list(ndeps = c(0.05,0.05)))
    MU <- est@coef[[1]]
    SI <- est@coef[[2]]
    
    ps <- vector()
    
    for (i in c(1:dim(bins)[1])) {
      ps <- c(ps, pnorm(bins[i,2], mean=MU, sd=SI, log.p=FALSE) - pnorm(bins[i,1], mean=MU, sd=SI, log.p=FALSE))
    }  

    #ps <- c(pnorm(bins[1,2], mean=MU, sd=SI, log.p=FALSE),
    #        pnorm(bins[2,2], mean=MU, sd=SI, log.p=FALSE) - pnorm(bins[2,1], mean=MU, sd=SI, log.p=FALSE),
    #        pnorm(bins[3,2], mean=MU, sd=SI, log.p=FALSE) - pnorm(bins[3,1], mean=MU, sd=SI, log.p=FALSE),
    #        pnorm(bins[4,2], mean=MU, sd=SI, log.p=FALSE) - pnorm(bins[4,1], mean=MU, sd=SI, log.p=FALSE),
    #        pnorm(bins[5,2], mean=MU, sd=SI, log.p=FALSE) - pnorm(bins[5,1], mean=MU, sd=SI, log.p=FALSE),
    #        pnorm(bins[6,1], mean=MU, sd=SI, log.p=FALSE, lower.tail=FALSE))
    
    total.count <- sum(bin.counts) / sum(ps)
  }

  ## Observation vector now has one entry for each bin, plus one entry for the estimated number of counts falling outside any of the bins
  o <- c(bin.counts, total.count - sum(bin.counts))

  est <- mle(minuslog=ll, 
           #start=list(mu=mu.i,std=sd.i),
           start=list(mu=MU,std=SI), 
           fixed=list(observations=o, bins=bins), 
           method="L-BFGS-B",
           lower=c(log10(minmean), minvar),
           # lower=c(0, 0), since we are in log space, I think there are no constraints on the mean and variance 
           upper=c(log10(maxmean), maxvar))#, 
           #control=list(ndeps = c(0.05,0.05)))
  print(length(est@coef))
  MU <- est@coef[[1]]
  SI <- est@coef[[2]]

  result <- c(MU, SI)
  names(result) <- c("mean", "sd")
  return(result)
}

## LOAD DATA 
loadReadCounts <- function(designfile, countsfile) {
  designDoc <- read.delim(designfile)
  counts <- read.delim(countsfile)
  
  # create a merged sheet
  mS <- merge(designDoc, counts)

  if (nrow(mS)==0) {
    writeLines("Unable to merge design doc and counts file", log)
    stop()
  }

  # if (! all(colnames(mS)[colnames(mS) %in% sort.params$bins$name] == sort.params$bins$name) ) {
  #   stop("Did not find columns matching bin names in the merged count table or did not find them in the same order")
  # }

  return(mS)
}

loadSortParams <- function(filename, mS, total.binname="Total") {
  ## Returns a list with the following items:
  ## bins: data.frame with the following columns:
  ##   name
  ##   mean (log10)
  ##   lowerBound (log10)
  ##   upperBound (log10)
  ##   count
  ## totalCount: Total count of cells sorted

  sort.params <- read.delim(filename)

  ## Check the sort params
  required.cols <- c("Mean","Bounds","Barcode","Count")
  if (! all(required.cols %in% colnames(sort.params)) ) {
    stop("Sort parameters file did not have the needed :", required.cols)
  }
  if (!(total.binname %in% sort.params$Barcode)) {
    stop(paste0("Barcode column in sort parameters file needs an entry called '",total.binname,"' to represent the total cell count from FACS"))
  }

  ## Extract info
  bin.indices <- which(sort.params$Barcode != total.binname)
  bin.names <- gsub("-",".",as.character(sort.params$Barcode[bin.indices]))
  bin.means <- log10(sort.params$Mean[bin.indices])
  bin.bounds <- do.call(rbind, lapply(strsplit(gsub("\\(| |\\)","", sort.params$Bounds[bin.indices]),","), as.numeric))
  total.count <- sort.params$Count[sort.params$Barcode == total.binname]
  seed <- log10(1+(sort.params[1,'Std.Dev.'] / sort.params[1,'Mean'])^2)

  if (! all(colnames(mS)[colnames(mS) %in% sort.params$bins$name] == sort.params$bins$name) ) {
    stop("Did not find columns matching bin names in the merged count table or did not find them in the same order")
  }
  
  filt.names <- vector()

  for (i in c(1:length(bin.names))) {
    name <- bin.names[i]
    if (sum(mS[,name]) > 0) {
      filt.names <- c(filt.names, bin.names[i])
    }
  } 

  bins <- data.frame(name=bin.names, mean=bin.means, lowerBound=log10(bin.bounds[,1]), upperBound=log10(bin.bounds[,2]), count=sort.params$Count[bin.indices], stringsAsFactors=F)
  rownames(bins) <- bin.names
  return(list(bins=bins[filt.names,], totalCount=total.count, sd.seed=seed))
}


rescaleReadCounts <- function(mS, sort.params, input="InputCount") {
  binNames <- sort.params$bins$name
  # write.table(mS[,binNames]/colSums(mS[,binNames]), file='/seq/lincRNA/Ben/base_editing/190415_BFF_spike_test/IL2RA-2/rs61839660-TATCTATTTTGGTCCCAAAC.NoSpike.BEFORE.txt', sep="\t", quote=F, row.names=F, col.names=T)
  mS[,binNames] <- sapply(binNames, function(binName) mS[,binName] / sum(mS[,binName]) * sort.params$bins[binName,"count"])
  # write.table(mS[,binNames], file='/seq/lincRNA/Ben/base_editing/190415_BFF_spike_test/IL2RA-2/rs61839660-TATCTATTTTGGTCCCAAAC.NoSpike.AFTER.txt', sep="\t", quote=F, row.names=F, col.names=T)
  mS[is.na(mS)] <- 0
  return(mS)
}


addWeightedAverage <- function(mS, sort.params) {
  mS$sum1 <- rowSums(mS[,sort.params$bins$name])
  mS$WeightedAvg <- ifelse(mS$sum1 == 0, 0, rowSums(t(sort.params$bins$mean * t(as.matrix(mS[,sort.params$bins$name])))) / mS$sum1)
  return(mS)
}

## Set up read count table and sort params
mS <- loadReadCounts(designDocLocation, countsLocation)
sort.params <- loadSortParams(sortParamsloc, mS)
mS <- rescaleReadCounts(mS, sort.params)
mS <- addWeightedAverage(mS, sort.params)
# write.table(mS, file='/seq/lincRNA/Ben/base_editing/190415_BFF_spike_test/IL2RA-2/rs61839660-TATCTATTTTGGTCCCAAAC.NoSpike.AFTER.txt', sep="\t", quote=F, row.names=F, col.names=T)

# run MLE
runMLE <- function(mS, sort.params) {
  ## seed mean and standard deviation
  sd.seed <- sort.params$sd.seed #.5 
  # print(sd.seed)
  bin.bounds <- as.matrix(sort.params$bins[,c("lowerBound","upperBound")])
  mleOuts <- data.frame(do.call(rbind, lapply(1:nrow(mS), function(i) 
    #tryCatch({ getNormalMLE(mS$WeightedAvg[i], sd.seed, mS[i,sort.params$bins$name], mS[i,inputCountCol], bin.bounds) }, error = function(err) { writeLines(paste0("MLE errored out at given initialization for guide: ", i, ", using weighted average instead"), log); result <- c(mS$WeightedAvg[i], sd.seed); names(result) <- c("mean", "sd"); return(result) }))))
    tryCatch({ getNormalMLE(mS$WeightedAvg[i], sd.seed, mS[i,sort.params$bins$name], mS[i,inputCountCol], bin.bounds) }, error = function(err) { write(paste0("MLE errored out at given initialization for guide: ", i, ", using weighted average instead: ", err),file=log,append=TRUE); result <- c(mS$WeightedAvg[i], sd.seed); names(result) <- c("mean", "sd"); return(result) }))))
  mS$logMean <- mleOuts$mean
  mS$logSD <- mleOuts$sd
  
  return(mS)
}

mS <- runMLE(mS, sort.params)

write.table(mS, file=outputmle, sep="\t", quote=F, row.names=F, col.names=T)
write("Test",file=log,append=TRUE)
#writeLines("Test", log)
