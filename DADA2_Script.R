#####Lake Okeechobee otros bac primers#####

#http://benjjneb.github.io/dada2/

setwd("~/Dropbox (UFL)/Lefler_Dissertation/LakeO_Phylo/")


######## PACKAGES & SUCH ###########
#BiocManager::install("ShortRead")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2")
#BiocManager::install("DECIPHER")
#BiocManager::install("Biostrings")
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#install.packages('ggthemes', dependencies = TRUE)
#install.packages('gg_ordiplot')
#install.packages("remotes")
#remotes::install_github("microbiome/microbiome")
#remotes::install_github("gauravsk/ranacapa")
#devtools::install_github("david-barnett/microViz")
#install.packages("remotes")
#remotes::install_github("jfq3/ggordiplots")
#BiocManager::install("ShortRead")
#BiocManager::install("Rhtslib")

#devtools::install_github("benjjneb/dada2")

library(readxl)
library('dada2')
library(microViz)
library(gclus)
library(ggrepel)
library(ggforce)
library(BiodiversityR)# also loads vegan
library(ranacapa)
library(microbiome)
library(tidyverse)
library(ggsci)
library(ggthemes)
library(ade4)
library(adegraphics)
library(vegan3d)
library(MASS)
library(ellipse)
library(FactoMineR)
library(rrcov)
library(missMDA)
library(Amelia)
library(BiodiversityR)
library(ggordiplots)
library(DECIPHER)
library(phyloseq)
library(Biostrings)
library(ShortRead)





###########################################################################
#####IMPORT READS


# CHANGE ME to the directory containing the fastq files after unzipping.
path <- "~/Dropbox (UFL)/Lefler_Dissertation/LakeO_Phylo/Reads"
list.files(path)

#####BEGIN HERE TO EXAMINE QUALITY OF READS

# Forward and reverse fastq filenames have format: C#_1.fastq and C#_2.fq.gz
fnFs <- sort(list.files(path, pattern="_1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Inspect read quality profiles
plotQualityProfile(fnFs[1:6]) #Forward Reads

plotQualityProfile(fnRs[1:6]) #Reverse Reads

#####BEGIN HERE FOR DADA2 PROCESS
###Filter and Trim

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


#Filter forward and reverse reads. TruncLen is determined by inspecting quality profiles and is dependent on sequences.
#In this example, forward reads are truncated at bp 285 and reverse reads are truncated at bp 275.
#maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE and maxEE=2 are the default parameters.
#The maxEE parameter sets the maximum number of "expected errors" allowed in a read, which is a better filter than simply averaging quality scores.

#Primers
#FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME # this is 515F
#REV <- "CCGYCAATTYMTTTRAGTTT"  ## CHANGE ME # this is 926YR


#the data looks great, therefore were are not going to trun anything
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, #truncLen=c(221,221),
                     trimLeft = c(19,20), #These are the lengths of the forward and reverse primers
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, verbose = TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)



#####Determining Error Rates#####
#Four options for learning error rates with NovaSeq data
#https://github.com/ErnakovichLab/dada2_ernakovichlab/tree/split_for_premise
#We went with option 4

#Option 1 from JacobRPrice alter loess arguments (weights and span and enforce monotonicity) benjjneb/dada2#1307
loessErrfun_mod1 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# check what this looks like
errF_1 <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod1,
  verbose = TRUE
)

errR_1 <- learnErrors(
  filtRs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod1,
  verbose = TRUE
)


#Option 2 enforce monotonicity only. Originally recommended in: benjjneb/dada2#791
loessErrfun_mod2 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        # mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}


# check what this looks like
errF_2 <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod2,
  verbose = TRUE
)


errR_2 <- learnErrors(
  filtRs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod2,
  verbose = TRUE
)

#Option 3 alter loess function (weights only) and enforce monotonicity From JacobRPrice benjjneb/dada2#1307
loessErrfun_mod3 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        # mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        # only change the weights
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot))
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# check what this looks like
errF_3 <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod3,
  verbose = TRUE
)



# check what this looks like
errR_3 <- learnErrors(
  filtRs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod3,
  verbose = TRUE
)

#Option 4 Alter loess function arguments (weights and span and degree, also enforce monotonicity) From Jonalim’s comment in benjjneb/dada2#1307

loessErrfun_mod4 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # jonalim's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),degree = 1, span = 0.95)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# check what this looks like
errF_4 <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod4,
  verbose = TRUE
)

errR_4 <- learnErrors(
  filtRs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod4,
  verbose = TRUE
)

#There are four options, none of which will yield “ideal” error plots. Instead look for the solution where the black line is continuously decreasing 
#(i.e. as quality scores improve on the x-axis the predicted error rate (y-axis) goes down) 
#and for plots that have points that mostly align with the black lines, although you will likely have some points along 0 on the y-axis.


# Original default recommended way (not optimal for NovaSeq data!)
errF_plot <- plotErrors(errF, nominalQ = TRUE)
errR_plot <- plotErrors(errR, nominalQ = TRUE)

errF_plot
errR_plot

# Trial 1 (alter span and weight in loess, enforce montonicity)
errF_plot1 <- plotErrors(errF_1, nominalQ = TRUE)
errR_plot1 <-plotErrors(errR_1, nominalQ = TRUE)

errF_plot1
errR_plot1

# Trial 2 (only enforce monotonicity - don't change the loess function)
errF_plot2 <- plotErrors(errF_2, nominalQ = TRUE)
errR_plot2 <-plotErrors(errR_2, nominalQ = TRUE)

errF_plot2
errR_plot2

# Trial 3 (alter loess (weights only) and enforce monotonicity)
errF_plot3 <- plotErrors(errF_3, nominalQ = TRUE)
errR_plot3 <-plotErrors(errR_3, nominalQ = TRUE)

errF_plot3
errR_plot3

# Trial 4 (alter loess (span, weight, and degree) and enforce monotonicity)
errF_plot4 <- plotErrors(errF_4, nominalQ = TRUE)
errR_plot4 <-plotErrors(errR_4, nominalQ = TRUE)

errF_plot4
errR_plot4



#this is how we used to do it,,,, dont do this lol
#The nbases value of 1e8 is for novaseq. 
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE) #Forward Reads
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE) #Reverse Reads

NSerrF_mon <- errF
NSerrR_mon <- errR


# assign any value lower than the Q40 probablity to be the Q40 value
#Forward Reads
NSnew_errF_out <- getErrors(NSerrF_mon) %>%
  data.frame() %>%
  mutate_all(funs(case_when(. < X40 ~ X40,
                            . >= X40 ~ .))) %>% as.matrix()
rownames(NSnew_errF_out) <- rownames(getErrors(NSerrF_mon))
colnames(NSnew_errF_out) <- colnames(getErrors(NSerrF_mon))
NSerrF_mon$err_out <- NSnew_errF_out
#Reverse Reads
NSnew_errR_out <- getErrors(NSerrR_mon) %>%
  data.frame() %>%
  mutate_all(funs(case_when(. < X40 ~ X40,
                            . >= X40 ~ .))) %>% as.matrix()
rownames(NSnew_errR_out) <- rownames(getErrors(NSerrR_mon))
colnames(NSnew_errR_out) <- colnames(getErrors(NSerrR_mon))
NSerrR_mon$err_out <- NSnew_errR_out





#####
#Applying sample inference algorithm
dadaFs <- dada(filtFs, err=errF_4, multithread=TRUE, pool=TRUE) #Forward Reads
dadaRs <- dada(filtRs, err=errR_4, multithread=TRUE, pool=TRUE) #Reverse Reads

dadaFs[[1]] #Identifies amount of sequence variants from forward reads.
dadaRs[[1]] #Identifies amount of sequence variants from reverse reads.

###Merge paired reads

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

###Construct sequence tables

seqtab <- makeSequenceTable(mergers) #Creates amplicon sequence table
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

##Optional (BUT NOT REALLY): remove non-target-length sequences from your sequence table

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 325:393] #get amplicons of the targeted length

###Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE) #combines a left-segment and a right-segment from two more abundant "parent" sequences.
dim(seqtab.nochim)

# Print percentage of our seqences that were not chimeric.
100*sum(seqtab.nochim)/sum(seqtab)

###Track reads through the pipeline. Determine the number of reads that made it through each step in the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
track

# tracking reads by percentage
track_pct <- track %>% 
  data.frame() %>%
  mutate(Sample = rownames(.),
         filtered_pct = ifelse(filtered == 0, 0, 100 * (filtered/input)),
         denoisedF_pct = ifelse(denoisedF == 0, 0, 100 * (denoisedF/filtered)),
         denoisedR_pct = ifelse(denoisedR == 0, 0, 100 * (denoisedR/filtered)),
         merged_pct = ifelse(merged == 0, 0, 100 * merged/((denoisedF + denoisedR)/2)),
         nonchim_pct = ifelse(nonchim == 0, 0, 100 * (nonchim/merged)),
         total_pct = ifelse(nonchim == 0, 0, 100 * nonchim/input)) 


# summary stats of tracked reads averaged across samples
track_pct_avg <- track_pct %>% summarize_at(vars(ends_with("_pct")), 
                                            list(avg = mean))
head(track_pct_avg)

track_pct_med <- track_pct %>% summarize_at(vars(ends_with("_pct")), 
                                            list(avg = stats::median))
head(track_pct_avg)
head(track_pct_med)

# Plotting each sample's reads through the pipeline
track_plot <- track %>% 
  data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  gather(key = "Step", value = "Reads", -Sample) %>%
  mutate(Step = factor(Step, 
                       levels = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim"))) %>%
  ggplot(aes(x = Step, y = Reads)) +
  geom_line(aes(group = Sample), alpha = 0.2) +
  geom_text(aes(label = Sample))+
  geom_point(alpha = 0.5, position = position_jitter(width = 0)) + 
  stat_summary(fun = median, geom = "line", group = 1, color = "steelblue", size = 1, alpha = 0.5) +
  stat_summary(fun = median, geom = "point", group = 1, color = "steelblue", size = 2, alpha = 0.5) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5), 
               geom = "ribbon", group = 1, fill = "steelblue", alpha = 0.2) +
  #geom_label(data = t(track_pct_avg[1:5]) %>% data.frame() %>% 
   #         rename(Percent = 1) %>%
    #         mutate(Step = c("filtered", "denoisedF", "denoisedR", "merged", "nonchim"),
     #              Percent = paste(round(Percent, 2), "%")),
      #     aes(label = Percent), y = 1.1 * max(track[,2])) +
  #geom_label(data = track_pct_avg[6] %>% data.frame() %>%
    #           rename(total = 1),
   #          aes(label = paste("Total\nRemaining:\n", round(track_pct_avg[1,6], 2), "%")), 
    #         y = mean(track[,6]), x = 6.5) +
  expand_limits(y = 1.1 * max(track[,2]), x = 7) +
  theme(legend.position = "none")+
  theme_classic()

track_plot

##Note: Outside of filtering, there should be no step in which a majority of reads are lost.
##If too many reads were lost return to filtering step.

#####END OF DADA2 PROCESSING
########################################################################


save.image(file="~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Jan2023/otros_physeq.RData") 





