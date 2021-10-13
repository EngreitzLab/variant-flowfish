diff --git a/workflow/rules/crispresso2.smk b/workflow/rules/crispresso2.smk
index fbc224f..0a48adc 100644
--- a/workflow/rules/crispresso2.smk
+++ b/workflow/rules/crispresso2.smk
@@ -7,14 +7,12 @@ rule run_crispresso:
 		read2=lambda wildcards: samplesheet.at[wildcards.SampleID,'fastqR2']
 	output:
 		directory('results/crispresso/CRISPResso_on_{SampleID}/')
-		#'crispresso/CRISPResso_on_{SampleID}/{AmpliconID}.Alleles_frequency_table_around_sgRNA_{GuideSpacer}.txt'
+		#'crispresso/CRISPResso_on_{SampleID}/{AmpliconID}.Alleles_frequency_table.txt'
 	params:
 		amplicon_id=lambda wildcards: samplesheet.at[wildcards.SampleID,'AmpliconID'],
 		amplicon_seq=lambda wildcards: samplesheet.at[wildcards.SampleID,'AmpliconSeq'],
-		guide=lambda wildcards: samplesheet.at[wildcards.SampleID,'GuideSpacer'],
 		q=config['crispresso_min_average_read_quality'],
 		s=config['crispresso_min_single_bp_quality'],
-		sample_cleavage=lambda wildcards: samplesheet.at[wildcards.SampleID, 'EditFromGuide']
 	#conda:
 	#    "envs/CRISPResso.yml"  
 	## 4/14/21 JE - Specifying the conda environment here is not working, and I am not sure why. Snakemake builds the conda environment, but then the conda environment doesn't work properly (CRISPResso not on the path)
@@ -33,11 +31,8 @@ rule run_crispresso:
 				--amplicon_seq {params.amplicon_seq} \
 				--amplicon_name {params.amplicon_id} \
 				--name {wildcards.SampleID} \
-				--guide_seq {params.guide} \
-				--cleavage_offset {params.sample_cleavage} \
-				--offset_around_cut_to_plot 7 \
 				-q {params.q} -s {params.s} || true'
 		"""
 
 # Ben's version also includes CRISPResso params:  --cleavage_offset -14 --offset_around_cut_to_plot 10 
-# To do: consider adding 				--max_paired_end_reads_overlap for FLASH overlap (calculate from read lengths, and length of amplicon)
\ No newline at end of file
+# To do: consider adding 				--max_paired_end_reads_overlap for FLASH overlap (calculate from read lengths, and length of amplicon)
diff --git a/workflow/rules/make_count_tables.smk b/workflow/rules/make_count_tables.smk
index 22acf7d..71d267e 100644
--- a/workflow/rules/make_count_tables.smk
+++ b/workflow/rules/make_count_tables.smk
@@ -5,7 +5,6 @@ import os
 import pandas as pd
 import numpy as np
 
-
 def make_count_table(samplesheet, group_col, group_id, bins, outfile, outfile_frequencies, samplesToExclude=None):
     ## Function to make a count table at various layers of resolution (e.g., by experiment, or by replicate, or by PCR replicate)
     ## To do: Move the python code for these rules into separate python scripts so they can be run independently of the snakemake pipeline (at least, this makes it easier to test and debug the code)
@@ -19,13 +18,11 @@ def make_count_table(samplesheet, group_col, group_id, bins, outfile, outfile_fr
     allele_tbls = []
 
     for idx, row in currSamples.iterrows():
-        file = "results/crispresso/CRISPResso_on_{SampleID}/{AmpliconID}.Alleles_frequency_table_around_sgRNA_{GuideSpacer}.txt".format(
-            SampleID=row['SampleID'], 
-            AmpliconID=row['AmpliconID'], 
-            GuideSpacer=row['GuideSpacer'])
+        file = "results/crispresso/CRISPResso_on_{SampleID}/Alleles_frequency_table.zip".format(
+            SampleID=row['SampleID'])
 
         if (os.path.exists(file)):
-            allele_tbl = pd.read_table(file)
+            allele_tbl = pd.read_csv(file)
             allele_tbl['#Reads'] = allele_tbl['#Reads'].astype(np.int32)
             ref_seq = allele_tbl.loc[allele_tbl['Aligned_Sequence'] == allele_tbl['Reference_Sequence'], 'Reference_Sequence'].values[0]
             allele_tbl = allele_tbl.loc[allele_tbl['Reference_Sequence'] == ref_seq] # necessary?
@@ -72,13 +69,11 @@ def make_flat_table(samplesheet, outfile):
 
     allele_tbls = []
     for idx, row in samplesheet.iterrows():
-        file = "{CRISPRessoDir}/{AmpliconID}.Alleles_frequency_table_around_sgRNA_{GuideSpacer}.txt".format(
-            CRISPRessoDir=row['CRISPRessoDir'],
-            AmpliconID=row['AmpliconID'],
-            GuideSpacer=row['GuideSpacer'])
+        file = "{CRISPRessoDir}/Alleles_frequency_table.zip".format(
+            CRISPRessoDir=row['CRISPRessoDir'])
 
         if (os.path.exists(file)):
-            allele_tbl = pd.read_table(file)
+            allele_tbl = pd.read_csv(file, sep='\t')
             allele_tbl['SampleID'] = row['SampleID']
             allele_tbls.append(allele_tbl)
 
diff --git a/workflow/scripts/AggregateDesiredAlleleCounts.R b/workflow/scripts/AggregateDesiredAlleleCounts.R
index 281142d..5aaa106 100644
--- a/workflow/scripts/AggregateDesiredAlleleCounts.R
+++ b/workflow/scripts/AggregateDesiredAlleleCounts.R
@@ -1,5 +1,5 @@
 # Jesse Engreitz
-# 5/30/21
+# 5/30/21 -- Katherine Guo updated 9/29/21
 # Rscript to aggregate allele frequency information across CRISPResso runs into a table for plotting and analysis
 
 
@@ -22,7 +22,6 @@ if (FALSE) {
   opt$variantCounts = "VariantCounts.flat.tsv.gz"
 }
 
-
 if (is.null(opt$outbase))
   stop("AggregateAlleleCounts: --outbase should be specified.\n")
 
@@ -39,17 +38,36 @@ if (opt$minFrequency < 0 | opt$minFrequency > 1)
 suppressPackageStartupMessages(library(tidyr))
 suppressPackageStartupMessages(library(dplyr))
 
+# Function to compare list of target sequences to match against a longer sequence. 
+# Return first match or 0 for no match.
+matchSeq <- function(seq, seqList) {
+  match = regmatches(seq, regexpr(paste(seqList, collapse="|"), seq)) %>% unlist()
+  # use gregexpr for all matches 
+  if (length(match) == 0) {
+    return("None") # HELP: if something doesn't match I think the sapply I use will break 
+  }
+  return(match)
+}
 
-getAlleleTable <- function(countsFlat, variantInfo, minFreq) {
-  countsFiltered <- merge(countsFlat, variantInfo %>% select(MappingSequence))
-
-  if (nrow(countsFiltered)==0) return(NULL)
 
+getAlleleTable <- function(countsFlat, variantInfo, minFreq) {
+  
+  # not sure if this is necessary, but does a similar function to the merge that was used previously
+  # countsFiltered <- filter(countsFlat, grepl(paste(variantInfo$MappingSequence, collapse="|"), MappingSequence)) 
+  # if (nrow(countsFiltered)==0) return(NULL)
+  countsFiltered <- countsFlat
+  
+  # apply matchSeq function to MappingSequence column and save matches in MatchSequence column
+  countsFiltered["MatchSequence"] <- sapply(countsFiltered["MappingSequence"], matchSeq, seqList=variantInfo$MappingSequence)
+  
   desiredCountsFlat <- countsFiltered %>% 
-    select(MappingSequence,`#Reads`,`%Reads`,SampleID) %>%         
-    group_by(MappingSequence,SampleID) %>% summarise(`%Reads`=sum(`%Reads`), `#Reads`=sum(`#Reads`)) %>% ungroup() %>% ## in some cases, CRISPResso reports the same Amplicon_Sequence twice ... if it aligns differently to reference
+    select(MatchSequence,`#Reads`,`%Reads`,SampleID) %>%         
+    group_by(MatchSequence,SampleID) %>% summarise(`%Reads`=sum(`%Reads`), `#Reads`=sum(`#Reads`)) %>% ungroup() %>% ## in some cases, CRISPResso reports the same Amplicon_Sequence twice ... if it aligns differently to reference
     as.data.frame()
 
+  desiredCountsFlat <- desiredCountsFlat %>%
+    rename(MappingSequence = MatchSequence) # rename MatchSequence back to MappingSequence so the merge with variantInfo at the end works 
+    
   desiredCountsTable <- desiredCountsFlat %>% select(-`#Reads`) %>% spread("SampleID","%Reads",fill=0) %>% as.data.frame()
 
   if (nrow(desiredCountsTable) > 1 & ncol(desiredCountsTable) > 2) {
@@ -61,7 +79,7 @@ getAlleleTable <- function(countsFlat, variantInfo, minFreq) {
     maxVal <- apply(desiredCountsTable[,-1], 1, max)
     desiredCountsTable <- desiredCountsTable[maxVal >= minFreq,]
   }
-
+  
   desiredCountsFlat <- merge(variantInfo %>% mutate(origOrder=1:n()), desiredCountsFlat) %>% arrange(origOrder) %>% select(-origOrder)
   desiredCountsTable <- merge(variantInfo %>% mutate(origOrder=1:n()), desiredCountsTable) %>% arrange(origOrder) %>% select(-origOrder)
 
@@ -69,9 +87,9 @@ getAlleleTable <- function(countsFlat, variantInfo, minFreq) {
 }
 
 
-countsFlat <- read.delim(gzfile(opt$variantCounts), check.names=F, stringsAsFactors=T) %>%
-            rename(MappingSequence=Aligned_Sequence) %>%
-            as.data.frame()
+countsFlat <- read.delim(gzfile(opt$variantCounts), check.names=F, stringsAsFactors=T)
+colnames(countsFlat)[1] <- "MappingSequence" 
+countsFlat <- as.data.frame(countsFlat)
 variantInfo <- read.delim(opt$variantInfo, check.names=F, stringsAsFactors=F)
 
 desiredCounts <- getAlleleTable(countsFlat, variantInfo, minFreq=opt$minFreq)
