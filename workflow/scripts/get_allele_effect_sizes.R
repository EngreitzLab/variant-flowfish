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
  # make_option(c("-dlo", "--designDocLocation"), type="character", help="Design document location, ex: ~/Documents/Harvard/LanderResearch/flo1/raw/160705.Ess.Design.txt"),
  make_option(c("-clo", "--countsLocation"), type="character", help="Counts document location, ex: ~/Documents/Harvard/LanderResearch/flo1/raw/4.tsv"),
  make_option(c("-sp", "--sortParamsloc"), type="character", help="MUST HAVE A MEAN, BOUNDS, AND BARCODE COLUMN!!, ex: ~/Documents/Harvard/LanderResearch/flo1/raw/sortParams/gata1ff3.txt"),
  make_option(c("-om", "--outputmle"), type="character", help="File to write the MLE mean into for each guide"),
  make_option(c("-rs", "--rescaledCounts"), type="character", help="File to check rescaled counts"),
  make_option(c("-l", "--log"), type="character", help="Log results of MLE"),
  make_option("--ignoreInputBinCounts", type="logical", default=FALSE, help="If TRUE, ignore the input bin bounds (Bin == 'All') in the MLE procedure")
)
opt <- parse_args(OptionParser(option_list=option.list))

## Set up variables
# designDocLocation <- opt$designDocLocation
countsLocation <- opt$countsLocation 
sortParamsloc <- opt$sortParamsloc
outputmle <- opt$outputmle
rescaledCounts <- opt$rescaledCounts
log <- opt$log 
useSeventhBin <- TRUE

print(paste("Counts Location: ", countsLocation))
print(paste("Sort Params Location: ", sortParamsloc))

MAXMEAN <- 10000 # in FACS space, max value a guide can take
MINMEAN <- 1  # in FACS space, max value a guide can take


## Save image with input parameters for debugging purposes
save.image(file=paste0(outputmle, ".RData"))


## load packages
# suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stats4))
suppressPackageStartupMessages(library(methods))
set.seed(100)

## FUNCTIONS

## One method of estimating the fraction of cells in the seventh bin
# takes as input the counts for allele i across the bins, the fraction of allele i in the input, and the total number of cells
# returns the number of "missing" cells (i.e. total * fraction - sum(observed))
estimateSeventhBinInput <- function(bin.counts, input.fraction, total.count) {
  return(total.count * input.fraction - sum(bin.counts))
}

## A different method of estimating the fraction of cells in the seventh bin
# takes as input a guess of starting mean+sd, bin counts and bin boundaries
# performs one step of MLE procedure, and uses this to estimate the fraction of cells that fall in 7th bin
# essentially performing one step of EM procedure
### TODO ###
# iterate EM until convergence
# katherine note: I added iteration to convergence but not totally sure this is right
estimateSeventhBinEM <- function(mu.guess, sd.guess, bin.counts, bins, minmean=MINMEAN, maxmean=MAXMEAN, minvar=0, maxvar=2) {
    o <- c(bin.counts, 0) # have to artificially pad with a 0 so the likelihood function doesn't complain
    prev_seventh_bin <- MAXMEAN
    seventh.bin.count <- 0
    iterations <- 0
    while((abs(seventh.bin.count - prev_seventh_bin) > 10) & (iterations < 100)) {
      iterations <- iterations + 1
      est <- mle(minuslog=ll,
              start=list(mu=mu.guess,std=sd.guess),
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
      prev_seventh_bin <- seventh.bin.count

      # use this to guess how many cells are in last bin
      seventh.bin.count <- sum(bin.counts) / sum(ps) * (1 - sum(ps))
      o <- c(bin.counts, seventh.bin.count)
    }

    return(seventh.bin.count)
}

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

  # initialize vector of probabilities of falling in bins 1..N
  pe <- vector()

  for (i in c(1:dim(bins)[1])) {
    pe <- c(pe, pnorm(bins[i,2], mean=mu, sd=std, log.p=FALSE) - pnorm(bins[i,1], mean=mu, sd=std, log.p=FALSE))
  }
  
  if (useSeventhBin){
    # add n+1th bin = p(fall outside bin)
    pe <- c(pe, 1-sum(pe)) # add a "bin" for the remaining cells
  }
  
  pe[pe <= 0] <- 10^-10 # remove any 0s from the probabilities (shouldn't happen but might)

  # assert that the lengths match
  if (length(pe) != length (observations)) {
    stop("Bin counts and bin boundaries don't have matching dimensions")
  }

  sum <- 0  # -log likelihood function
  for (i in c(1:length(observations))) {
    sum <- sum - (log(pe[i]) * observations[i]) # log(p^k)
  }

  return(sum)
}

# Wrap the maximum likelihood estimator
# getNormalMLE <- function(mu.i, sd.i, bin.counts, total.count, bins, minmean=MINMEAN, maxmean=MAXMEAN, minvar=0, maxvar=1) {
# getNormalMLE <- function(mu.i, sd.i, bin.counts, seventh.bin.count, bins, minmean=MINMEAN, maxmean=MAXMEAN, minvar=0, maxvar=1) {
getNormalMLE <- function(mu.i, sd.i, bin.counts, bins, input.present, idx, total.count, mS, minmean=MINMEAN, maxmean=MAXMEAN, minvar=0, maxvar=2) {
  ## mu.i = initial guess of the mean
  ## si.i = initial guess of the standard deviation
  ## bin.counts = vector of counts per bin (observed)
  ## bins = matrix of bin boundaries (nrows = nbins, ncols = 2)
  ## returns mean and standard deviation in log-space
  if (idx %% 1000 == 0) cat(paste0("get_allele_effect_sizes.R: Processed ", idx, " alleles.\n"))

  # cast to numeric vector
  bin.counts <- as.numeric(as.matrix(bin.counts))

  if (useSeventhBin) {
    # first, need to figure out how to estimate the seventh bin counts
    seventh.bin.count <- 0
    if (input.present) {
      seventh.bin.count <- estimateSeventhBinInput(bin.counts, mS$input.fraction[idx], total.count)
    }

    if (seventh.bin.count <= 0 || !input.present) {
      print("Estimating seventh bin count with MLE")
      seventh.bin.count <- estimateSeventhBinEM(mu.i, sd.i, bin.counts, bins)
    }

    # add on "seventh" bin counts
    # o <- c(bin.counts, total.count - sum(bin.counts))
    o <- c(bin.counts, seventh.bin.count)
  } else {
    o <- bin.counts
  }
  ## Observation vector now has one entry for each bin, plus one entry for the estimated number of counts falling outside any of the bins (n+1 = "7")
  ## bins still has only n (=6) entries, for the observed bins, but the log-likelihood function `ll` adds in the n+1th bin (7th bin = 1-sum(p_i))

  est <- mle(minuslog=ll,
           start=list(mu=mu.i,std=sd.i),
           # start=list(mu=MU,std=SI),
           fixed=list(observations=o, bins=bins),
           method="L-BFGS-B",
           lower=c(log10(minmean), minvar),
           # lower=c(0, 0), since we are in log space, I think there are no constraints on the mean and variance
           upper=c(log10(maxmean), maxvar))#,
           #control=list(ndeps = c(0.05,0.05)))
  MU <- est@coef[[1]]
  SI <- est@coef[[2]]

  result <- c(MU, SI)
  names(result) <- c("mean", "sd")
  return(result)
}

## LOAD DATA
# loadReadCounts <- function(designfile, countsfile) {
loadReadCounts <- function(countsfile) {
  # designDoc <- read.delim(designfile)
  counts <- read.delim(countsfile)

  # create a merged sheet
  # mS <- merge(designDoc, counts)

  if (nrow(counts)==0) {
    writeLines("Unable to merge design doc and counts file", log)
    stop()
  }

  # if (! all(colnames(mS)[colnames(mS) %in% sort.params$bins$name] == sort.params$bins$name) ) {
  #   stop("Did not find columns matching bin names in the merged count table or did not find them in the same order")
  # }

  return(counts)
}


getBinNames <- function(counts) {
  col.names <- colnames(counts)
  bin.names <- col.names [! col.names %in% c('AmpliconID', 'MappingSequence','All', 'Neg', 'Reference_Name', 'Match_Sequence', 'VariantID', 'RefAllele')]
  return(bin.names)
}



loadSortParams <- function(filename, counts, binNames, totalBinName="Total") {
  ## Function to read in sort params files, including detecting among
  ## several possible formats from different sorting instruments
  ##
  ## To do: Consider moving this into a separate script
  ## To do: Consolidate so that we need either counts or binNames but not both.  

  sort.params <- read.csv(filename, sep=',')

  print("Detecting sort params file format...")

  required.cols.astrios <- c("Mean","Bounds","Barcode","Count")
  required.cols.bigfoot <- c("Mean","Min","Max","Barcode","Count","StdDev")
  required.cols.influx <- c("Mean","Min","Max","Barcode","Count")

  if (all(required.cols.astrios %in% colnames(sort.params))) {
    print("  Attempting to load sort parameters using Astrios file format...")
    result <- loadSortParams_Astrios(filename, counts, total.binname=totalBinName)
    print("  Sort params:")
    print(result)
  } else if (all(required.cols.bigfoot %in% colnames(sort.params))) {
    print("  Attempting to load sort parameters using BigFoot file format...")
    result <- loadSortParams_BigFoot(filename, binNames, total.binname=totalBinName)
    print("  Sort params:")
    print(result)
  } else if (all(required.cols.influx %in% colnames(sort.params))) {
    print("  Attempting to load sort parameters using Influx file format...")
    result <- loadSortParams_Influx(filename, binNames, total.binname=totalBinName)
    print("  Sort params:")
    print(result)
  } else {
    stop("ERROR: Did not find required columns for any of the sort params file formats.\n")
  }

  return(result)
}


loadSortParams_Astrios <- function(filename, mS, total.binname="Total") {
  ## Returns a list with the following items:
  ## bins: data.frame with the following columns:
  ##   name
  ##   mean (log10)
  ##   lowerBound (log10)
  ##   upperBound (log10)
  ##   count
  ## totalCount: Total count of cells sorted

  sort.params <- read.csv(filename, strip.white=TRUE)

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
  if (length(filt.names) == 0) stop("Stopping because no bins have counts in the input count table.\n")

  bins <- data.frame(name=bin.names, mean=bin.means, lowerBound=log10(bin.bounds[,1]), upperBound=log10(bin.bounds[,2]), count=sort.params$Count[bin.indices], stringsAsFactors=F)
  rownames(bins) <- bin.names
  return(list(bins=bins[filt.names,], totalCount=total.count, sd.seed=seed))
}



loadSortParams_BigFoot <- function(filename, bin.names, total.binname="Total", full.file=TRUE) {
  ## Returns a list with the following items:
  ## bins: data.frame with the following columns:
  ##   name
  ##   mean (log10)
  ##   lowerBound (log10)
  ##   upperBound (log10)
  ##   count
  ## totalCount: Total count of cells sorted

  sort.params <- read.csv(filename, strip.white=TRUE)

  ## Check the sort params
  # check that it has all the required columns
  required.cols <- c("Mean","Min","Max","Barcode","Count","StdDev")
  if (! all(required.cols %in% colnames(sort.params)) ) {
    stop("Sort parameters file ", filename, " did not have the needed required column(s): ", setdiff(required.cols,colnames(sort.params)))
  }
  # check that it has a total count
  if (!(total.binname %in% sort.params$Barcode)) {
    stop(paste0("Barcode column in sort parameters file ", filename, " needs an entry called '",total.binname,"' to represent the total cell count from FACS"))
  }
  # check that it has all the bins that the countsFile has
  if (sum(sort.params$Barcode %in% bin.names) != length(bin.names)) {
    print(bin.names)
    print(sort.params$Barcode)
    stop("Sort parameters file ", filename, " did not have all the bins")
  }

  ## Extract info from sortParams for each bin
  bin.indices <- which(sort.params$Barcode %in% bin.names)
  bin.names.inorder <- sort.params$Barcode[bin.indices]
  # bin.names.inorder <- gsub("-",".",as.character(sort.params$Barcode[bin.indices]))

  if (full.file) {
    # need to treat the rows differently and cast them away from "Factors"
    bin.means <- log10(as.numeric(as.character(sort.params$Mean)[bin.indices]))  ## FOR BIGFOOT, this number appears to be wrong
    bin.mins <- log10(as.numeric(as.character(sort.params$Min)[bin.indices]))
    bin.maxs <- log10(as.numeric(as.character(sort.params$Max)[bin.indices]))
    seed <- log10(1+(as.numeric(as.character(sort.params$StdDev[sort.params$Barcode == total.binname])) / as.numeric(as.character(sort.params$Mean[sort.params$Barcode == total.binname])))^2)
  } else {
    bin.means <- log10(sort.params$Mean[bin.indices])
    bin.mins <- log10(sort.params$Min[bin.indices])
    bin.maxs <- log10(sort.params$Max[bin.indices])
    seed <- log10(1+(sort.params$StdDev[sort.params$Barcode == total.binname] / sort.params$Mean[sort.params$Barcode == total.binname])^2)
  }

  bin.counts <- sort.params$Count[bin.indices]

  mu.seed <- median(bin.mins)

  # compute some extra numbers which will help our estimation
  total.count <- sort.params$Count[sort.params$Barcode == total.binname]

  bins <- data.frame(name=bin.names.inorder, mean=bin.means, lowerBound=bin.mins, upperBound=bin.maxs, count=bin.counts, stringsAsFactors=F)
  rownames(bins) <- bin.names.inorder
  return(list(bins=bins, totalCount=total.count, sd.seed=seed, mu.seed=mu.seed))
}


loadSortParams_Influx <- function(filename, bin.names, total.binname="Total", full.file=TRUE) {
  ## 7/26/21 temporary fix for reading Influx data
  ## TODO: Factor out these different loadSortParams files and add code to fix these before launching the pipeline

  ## Returns a list with the following items:
  ## bins: data.frame with the following columns:
  ##   name
  ##   mean (log10)
  ##   lowerBound (log10)
  ##   upperBound (log10)
  ##   count
  ## totalCount: Total count of cells sorted

  sort.params <- read.csv(filename, strip.white=TRUE)

  ## Check the sort params
  # check that it has all the required columns
  required.cols <- c("Mean","Min","Max","Barcode","Count")
  if (! all(required.cols %in% colnames(sort.params)) ) {
    stop("Sort parameters file ", filename, " did not have the needed required column(s): ", setdiff(required.cols,colnames(sort.params)))
  }
  # check that it has a total count
  if (!(total.binname %in% sort.params$Barcode)) {
    stop(paste0("Barcode column in sort parameters file ", filename, " needs an entry called '",total.binname,"' to represent the total cell count from FACS"))
  }
  # check that it has all the bins that the countsFile has
  if (sum(sort.params$Barcode %in% bin.names) != length(bin.names)) {
    print(bin.names)
    print(sort.params$Barcode)
    stop("Sort parameters file ", filename, " did not have all the bins")
  }

  ## Extract info from sortParams for each bin
  bin.indices <- which(sort.params$Barcode %in% bin.names)
  bin.names.inorder <- sort.params$Barcode[bin.indices]
  # bin.names.inorder <- gsub("-",".",as.character(sort.params$Barcode[bin.indices]))

  if (full.file) {
    # need to treat the rows differently and cast them away from "Factors"
    bin.means <- log10(as.numeric(as.character(sort.params$Mean)[bin.indices]))
    bin.mins <- log10(as.numeric(as.character(sort.params$Min)[bin.indices]))
    bin.maxs <- log10(as.numeric(as.character(sort.params$Max)[bin.indices]))
    seed <- log10(1+median(as.numeric(as.character(sort.params$Mean)[bin.indices])))-1 ## For influx, missing this data — make a guess
  } else {
    bin.means <- log10(sort.params$Mean[bin.indices])
    bin.mins <- log10(sort.params$Min[bin.indices])
    bin.maxs <- log10(sort.params$Max[bin.indices])
    seed <- log10(1+median(as.numeric(as.character(sort.params$Mean)[bin.indices])))-1 ## For influx, imssing this data - make a guess
  }

  bin.counts <- sort.params$Count[bin.indices]

  mu.seed <- median(bin.mins)

  # compute some extra numbers which will help our estimation
  total.count <- sort.params$Count[sort.params$Barcode == total.binname]

  bins <- data.frame(name=bin.names.inorder, mean=bin.means, lowerBound=bin.mins, upperBound=bin.maxs, count=bin.counts, stringsAsFactors=F)
  rownames(bins) <- bin.names.inorder
  return(list(bins=bins, totalCount=total.count, sd.seed=seed, mu.seed=mu.seed))
}




rescaleReadCounts <- function(mS, sort.params, bin.names, input="InputCount") {
  binNames <- bin.names
  # write.table(mS[,binNames]/colSums(mS[,binNames]), file='/seq/lincRNA/Ben/base_editing/190415_BFF_spike_test/IL2RA-2/rs61839660-TATCTATTTTGGTCCCAAAC.NoSpike.BEFORE.txt', sep="\t", quote=F, row.names=F, col.names=T)
  mS[,binNames] <- sapply(binNames, function(binName) mS[,binName] / sum(mS[,binName]) * sort.params$bins[binName,"count"])
  # write.table(mS[,binNames], file='/seq/lincRNA/Ben/base_editing/190415_BFF_spike_test/IL2RA-2/rs61839660-TATCTATTTTGGTCCCAAAC.NoSpike.AFTER.txt', sep="\t", quote=F, row.names=F, col.names=T)
  mS[is.na(mS)] <- 0
  return(mS)
}


addWeightedAverage <- function(mS, sort.params, bin.names) {
  mS$sum1 <- rowSums(mS[,bin.names])
  mS$WeightedAvg <- ifelse(mS$sum1 == 0, 0, rowSums(t(sort.params$bins[bin.names,'mean'] * t(as.matrix(mS[,bin.names])))) / mS$sum1)
  return(mS)
}

## Set up read count table and sort params
counts <- loadReadCounts(countsLocation)
bin.names <- getBinNames(counts)
sort.params <- loadSortParams(sortParamsloc, counts, bin.names)
counts <- rescaleReadCounts(counts, sort.params, bin.names)
counts <- addWeightedAverage(counts, sort.params, bin.names)

# run MLE
runMLE <- function(mS, sort.params, ignoreInputBinCounts=FALSE) {
  ## seed mean and standard deviation
  mu.seed <- sort.params$mu.seed
  sd.seed <- sort.params$sd.seed #.5
  # print(sd.seed)
  bin.bounds <- as.matrix(sort.params$bins[,c("lowerBound","upperBound")])
  total.count <- sort.params$totalCount

  input.bin.name <- 'All'
  input.present <- FALSE

  if (input.bin.name %in% colnames(mS) & !ignoreInputBinCounts) {
    mS$input.fraction <- mS[,input.bin.name] / sum(mS[,input.bin.name])
    input.present <- TRUE
  }

  ############# MAKE SURE BIN.BOUNDS AND mS[i,bin.names] ARE SORTED IN THE SAME ORDER!!! #############
  mleOuts <- data.frame(do.call(rbind, lapply(1:nrow(mS), function(i)
    #tryCatch({ getNormalMLE(mS$WeightedAvg[i], sd.seed, mS[i,sort.params$bins$name], mS[i,inputCountCol], bin.bounds) }, error = function(err) { writeLines(paste0("MLE errored out at given initialization for guide: ", i, ", using weighted average instead"), log); result <- c(mS$WeightedAvg[i], sd.seed); names(result) <- c("mean", "sd"); return(result) }))))
    # tryCatch({ getNormalMLE(mS$WeightedAvg[i], sd.seed, mS[i,sort.params$bins$name], mS[i,inputCountCol], bin.bounds) }, error = function(err) { write(paste0("MLE errored out at given initialization for guide: ", i, ", using weighted average instead: ", err),file=log,append=TRUE); result <- c(mS$WeightedAvg[i], sd.seed); names(result) <- c("mean", "sd"); return(result) }))))
    # estimate the number of cells falling in n+1th bin
    # estimated.seventh.bin <- estimateSeventhBinInput(mS[i,bin.names], mS[i,input.bin.name] / sum(mS[,input.bin.name]), total.count)
    tryCatch({
      getNormalMLE(mu.seed, sd.seed, mS[i,bin.names], bin.bounds, input.present, i, total.count, mS)
    }, error = function(err) {
      write(paste0("MLE errored out at given initialization for guide: ", i, ", using weighted average instead: ", err),file=log,append=TRUE)
      result <- c(mS$WeightedAvg[i], sd.seed); names(result) <- c("mean", "sd")
      return(result)
    }))))

  mS$logMean <- mleOuts$mean
  mS$logSD <- mleOuts$sd

  return(mS)
}


# write(pe[7],file=log,append=TRUE)
mleOutput <- runMLE(counts, sort.params, opt$ignoreInputBinCounts)

write.table(mleOutput, file=outputmle, sep="\t", quote=F, row.names=F, col.names=T)
write.table(counts, file=rescaledCounts, sep="\t", quote=F, row.names=F, col.names=T)
write("Finished",file=log,append=TRUE)
