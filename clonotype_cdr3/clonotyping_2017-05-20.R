# 2017-04-12:   Added output summary file <yyyy-mm-dd>_test_summary.txt 
#               by sink()'ing every cat() expression.
#               Outputs are named with sample name and date now
# future:
# put everything in a big function with subfunctions, assign path names
# include 2nd function with graphing capabilities
# first step function, merge dataframes 1_Summary and 5_AA

wkdir("./clonotype_cdr3/test/")
rm(list=ls())
library(stringdist)
library(compiler)
library(tidyverse)
library(Biostrings)
require(plyr) 


#' NOTES: Install dendextend to fix the code.
#' Original problem: The same sequences were grouped into different clusters
#' Fix: Use Cutree from dextend package (also improved clustering algorithm )
#' Added spike-in filter for human and mouse

if(!require(dendextend)){
  installe.packages("dendextend")
  library(dendextend)
}
if(!require(seqinr)){
    install.packages('seqinr')
    require(seqinr)
}

## Read spike list
fpath_spikes <- paste("/Users/mfilipav/seqdata/spike_ins/spike_cdr3_aa.csv")
spikes <- read.csv(fpath_spikes)
######################################################################################################
# reads IMGT files
sample_name <- "bcma7"
filename_prefix <- paste0(Sys.Date(),"_", sample_name, "_")

setwd("./clonotype_cdr3/test/")
# create summary file:
sink(paste0(filename_prefix, "summary.txt"))

#writeLines(filename_prefix, summary_file, sep="\n")
#fpath <- paste("~/seqdata/raw_data/aligned_april2017/20170413_Analysis/bcma_6/Outfiles/IMGT_vdj/5_AA-Outfiles.txt, sep ="")

#fpath2 <- paste("1_Summary_.imgt.txt", sep ="")
imgt_5aa_original <- read.delim("./clonotype_cdr3/test/5_AA-Outfiles.txt")
imgt_1sum <- read.delim("./clonotype_cdr3/test/1_Summary_Outfiles.txt")
imgt_1sum <- imgt_1sum[,c(1, 31)]

initial_reads <- length(imgt_1sum[,1])
cat(paste("Initial reads into analysis:", initial_reads), "\n")
#writeLines(paste("Initial reads into analysis:", initial_reads), summary_file, sep="\n")
?writeLines
# Merge two IMGT files, write output
imgt_5aa<- merge(imgt_5aa_original, imgt_1sum)

# 1. DATA CLEAN-UP
# Filter for productive V-genes
imgt_5aa_prod <- imgt_5aa[imgt_5aa$Functionality == "productive",]
productive_reads <- length(imgt_5aa_prod[,1])
cat(paste("Reads after removing Unproductive sequences:", productive_reads), "\n")
#write(paste("Reads after removing Unproductive sequences:", productive_reads), summary_file, append=TRUE)


## Filter for at least two occurances for each aa VDJ seq
imgt_5aa_prod_filtered <-subset(imgt_5aa_prod, duplicated(imgt_5aa_prod$V.D.J.REGION))
# Remove empty VDJ and CDR3 column values
imgt_5aa_prod_filtered <- subset(
    imgt_5aa_prod_filtered, subset = imgt_5aa_prod_filtered$V.D.J.REGION != ""
    )
imgt_5aa_prod_filtered <- subset(
    imgt_5aa_prod_filtered, subset = imgt_5aa_prod_filtered$CDR3.IMGT != ""
    )

productive_reads_flitered <- length(imgt_5aa_prod_filtered[,1])
cat(paste("Reads for which VDJ aa sequence appears at least twice:", productive_reads_flitered), "\n")



# Write out the cleaned up data
write.table(imgt_5aa_prod_filtered, file = paste0(filename_prefix, "1S_5AA_merged_cleaned.txt"), sep = "\t")


# 2. INITIAL STATS
# total number of different VDJ genes
# do VDJ rank table now
# number of different VDJ genes
# number of different VDJ genes above threshold (5)
# V-gene allele frequency plot
# CDR3 lengths plot
####

# # v-gene frequencies
# vdj_frequencies <- sort(prop.table(table(imgt_5aa_prod_filtered$V.GENE.and.allele))*100, decreasing=T)
# class(vdj_frequencies)
# vdj_frequencies_significant <- as.data.frame(vdj_frequencies)
# vdj_frequencies_significant <- vdj_frequencies_significant[which(vdj_frequencies_significant$Freq > .1), ]
# print("VDJ gene frequencies:")
# vdj_frequencies_significant
# #vdj_frequencies_significant<- as.table(vdj_frequencies_significant$Var1, vdj_frequencies_significant$Freq)
# #class(vdj_frequencies_significant)
# #barplot(vdj_frequencies_significant, ylab = "Frequency (%)")
# 
# 
# # total number of different VDJ genes
# vdj_threshold <- 100
# vdj_table <- table(imgt_5aa_prod_filtered$V.D.J.REGION)
# 
# df <- as.data.frame(vdj_table)
# names(df) <- c("VDJ","count") #maybe try to fix name here?!
# df$rank <- rank(-df$count, na.last = NA, ties.method = "min")
# vdjrank <- df[order(df$rank,decreasing = F),]
# vdjrank <- vdjrank[ which(vdjrank$count > 0), ] # trim all seqs with less than 5 CDR3 aa seqs or 5 VH aa seqs
# View(vdjrank)
# unique_vdj <- nrow(vdjrank)
# unique_vdj_threshold <- nrow(vdjrank[ which(vdjrank$count > vdj_threshold), ])
# cat(paste("Unique VDJ aa seqs:", unique_vdj))
# cat(paste("Unique VDJ aa seqs that have more than", vdj_threshold, "counts:", unique_vdj_threshold))
# 
# vdjrank_top50 <- vdjrank[ which(vdjrank$rank < 51), ]
# plot(x = vdjrank_top50$rank, y = vdjrank_top50$count, xlab = "rank", ylab = "VDJ aa counts")
# write.table(vdjrank, file = "vdjranks_26032017.txt", sep = "\t")
# 


# 3. CDR3 AA RANKING
## Ranking unique CDR3 aa by abundance
cdr3_threshold <- 0

cdr3table <- table(imgt_5aa_prod_filtered$CDR3.IMGT)

df <- as.data.frame(cdr3table)
names(df) <- c("CDR3","count") #maybe try to fix name here?!
df$rank <- rank(-df$count, na.last = NA, ties.method = "min")
cdr3rank <- df[order(df$rank,decreasing = F),]
cdr3rank <- cdr3rank[ which(cdr3rank$count > 0), ] # trim all seqs with less than 5 CDR3 aa seqs or 5 VH aa seqs

unique_cdr3 <- nrow(cdr3rank)
unique_cdr3_threshold <- nrow(cdr3rank[ which(cdr3rank$count > cdr3_threshold), ])

cat(paste("Unique HCDR3 aa seqs:", unique_cdr3) , "\n")
cat(paste("Unique HCRD3 aa seqs that have more than", cdr3_threshold, "counts:", unique_cdr3_threshold), "\n")

cdr3rank_top50 <- cdr3rank[ which(cdr3rank$rank < 51), ]
plot(x = cdr3rank_top50$rank, y = cdr3rank_top50$count, xlab = "rank", ylab = "HCDR3 aa counts")

write.table(cdr3rank, file = paste0(filename_prefix,"cdr3_ranking_threshold_", 
                                    cdr3_threshold, ".txt"),
            sep = "\t")

#spec accum plot
unique_appearance_address <- seq_along(imgt_5aa_prod_filtered$CDR3.IMGT)[cdr3rank$CDR3]
sort(unique_appearance_address)
plot(x = unique_appearance_address)



## Clonotyping by % seq identity
## In: unique cdr3 aa seqs ranked by abundance
## Out: clonotype clusters dataframe with total clonotype count number
clone_count_threshold <-1
dataframe <- cdr3rank[ which(cdr3rank$count > clone_count_threshold), ]
cdr3s <- dataframe$CDR3
counts <- dataframe$count

#' seq = CDR3s, #labels = V/J for example, or any other qualify by which you 
#' would like to divide CDR3s prior to clonotyping, 
#' #similarity = sequence similarity 
#' e.g., 0.1 means clonotype by 90% similarity)
clonotyping_function <- function(seq, labels, similarity){ 
  
  tapply(seq, labels, function (x) {
    cat("c")
    
    if(length(x)>1){ 
      maxL <- function(x) {
        i <- unlist(lapply(c(2:length(x)), function(y) seq(y, length(x), 1)))
        j <- rep(1:(length(x)-1), rev(1:(length(x)-1)))
        pmax(nchar(x[i]), nchar(x[j]))
      }

      ## faster way to calculate distance clustering (but not the fastest )
      clust <- hclust(as.dist(stringdistmatrix(x,  method = 'lv'))/maxL(x))
      clust$labels <- x
      
      clust <- as.dendrogram(clust)
      
      ## cut the cluster dendrogram by similarity
      cut_cdr3s <- cutree(clust, h = similarity) 
      
      ## group clones in list by similarity
      tapply(names(cut_cdr3s), cut_cdr3s,  function(x)x) 
    }
    else{
      x
    } 
  })
}

## makes function a bit faster/more efficient
clonotyping_function <- cmpfun(clonotyping_function) 

#' Performs clonotyping by 75% similarity. 
#' Since you wished to disregard VJ information,
#' I gave all CDR3s the arbitraty label "1???.
clonotyping_threshold <- 0.25
cdr3_clonotyping_list <- clonotyping_function(as.character(cdr3s), rep("1", length(cdr3s)), 
                                              clonotyping_threshold)  

# Counts the number of clones per cluster
count_list <- vector("list", length(cdr3_clonotyping_list[[1]]) )   
for(i in 1:length(cdr3_clonotyping_list[[1]])) {
    
    ## catch the clones from first clonotype, [[1]][[2]] for 2nd clonotype
    test <- c(cdr3_clonotyping_list[[1]][[i]])      
    
    ## gives sum of clone counts from one clonotype
    cluster_sum <- sum(dataframe[cdr3s %in% test,"count"])  
    
    ## add clonotype counts to list,
    ## update list index, until number = # clonotypes
    count_list[i]  <- cluster_sum
}

# Transforms data into a df
l <- lapply(cdr3_clonotyping_list[[1]], as.data.frame)
final_data_df <- as.data.frame(rbind.fill(lapply(l, function(x)as.data.frame(t(x)))))

final_data_df <- data.frame(lapply(final_data_df, as.character), stringsAsFactors=FALSE)

clone_counts <- unlist(count_list)

final_df <- cbind(clone_counts, final_data_df)
number_clonotypes <- nrow(final_df) 
cat(paste("Number of clonotypes at h =", clonotyping_threshold, "aa sequence identity is:", number_clonotypes, "\n"))
cat(paste("Clones must appear more than", clone_count_threshold, "times to be considered for clonotyping", "\n"))
    
sum(final_df$clone_counts)
######
# clean up human and mouse spike-ins using a spike_cdr3_aa.csv

#mask <- paste('AR', spikes_hum$`spike_cdr3_aa$Sequence`, sep ='') %in% final_df$V1
mask <- final_df$V1%in%paste('AR', spikes$Sequence, sep ='')
final_df <- final_df[!mask, ]
sum(final_df$clone_counts)

View(final_df)
write.table(final_df, file = paste0(filename_prefix,"cdr3_clonotypes_h_", 
                                    clonotyping_threshold, ".txt"), sep = "\t")



sink()
# SINK ENDS HERE
################################################################################
# 3. Plot clonotype CDR3 aa lengths
cdr3_clonotypes <- as.vector(final_df[,2])
cdr3_clonotype_lengths <- sapply(cdr3_clonotypes, nchar)
class(cdr3_clonotypes)

par(las=2) # make label text perpendicular to axis
par(mar=c(15,4,4,6), xpd = TRUE)
plot(cdr3_clonotype_lengths, type="o", col="dark blue", lwd=2, 
    xlab=".", ylab="CDR3 length (a.a.)", ylim=c(0,max(nchar(cdr3_clonotypes))),
    main="Length of ROR1 HCDR3 clonotype aa sequences")
axis(1, at=NULL, labels = FALSE )
text(seq(1, length(cdr3_clonotypes), by=1), par("usr")[1] - 3, 
     labels = cdr3_clonotypes, srt = 90, pos = 2, xpd = TRUE)
main="Length of unique ROR1 HCDR3 aa sequences used in clonotyping"


## Clonotypes CDR3 aa constitution color box plots
seqs <- cdr3_clonotypes
length(seqs)
AAs <- alphabetFrequency(AAStringSet(seqs)) ## alphabetFrequency also works on sets
AAs <- AAs[,1:20]
layout(matrix(c(1,2), 1, 2, byrow = TRUE), width=c(0.7,0.3)) ## allow to place the legend
nAAs = AAs/rowSums(AAs) ## if you want to normalize by total length, use this
par(las=2) # make label text perpendicular to axis
par(mar=c(13,4,4,5)) # increase y-axis margin.
barplot(t(nAAs),col=rainbow(ncol(AAs)), main = "CDR3 aa composition of the dominant clone from each clonotype", 
        names.arg = seqs, bty = "L" ) ## yields a stacked barplot, one bar per sequence
par(xpd=TRUE)
legend("topright", inset=c(-0.2, 0), legend=rev(colnames(AAs)), fill=rev(rainbow(ncol(AAs))),
       cex=0.8, xjust = 2)



## Unique clone CDR3 aa lengths
cdr3_clones <- cdr3rank[ which(cdr3rank$count > 0), ]
cdr3_clones <- as.character(cdr3_clones$CDR3, mode = "list")    # into character
cdr3_clone_lengths <- sapply(cdr3_clones, nchar)
plot(cdr3_clone_lengths, xlab="Clone rank (by abundance)",
     ylab="CDR3 length (a.a.)", col="dark blue", 
     lwd=2, ylim=c(0,max(nchar(cdr3_clones)))) # add type="o" for lines
title(main="Length of unique ROR1 HCDR3 aa sequences used in clonotyping", 
      col.main="black", font.main=1)

## Density plot of unique HCDR3 aa clone length distribution
plot(density(cdr3_clone_lengths), xlab="CDR3 length (a.a)", main="", 
     col="dark blue", lwd=5)
title(main="Density distribution of unique ROR1 HCDR3 aa sequence lengths", 
      col.main="black", font.main=1)

## Density plot of all quality filtered HCDR3 aa clones' lengths distro
cdr3_seqs <- as.character(imgt_5aa_prod_filtered$CDR3.IMGT)
summary(cdr3_seqs)
cdr3_seqs_lengths <- sapply(cdr3_seqs, nchar)
plot(density(cdr3_seqs_lengths), xlab="CDR3 length (a.a)", main="", 
     col="dark blue", lwd=3)
title(main="Density distribution of HCDR3 aa sequence lengths from all ROR1 reads", 
      col.main="black", font.main=1)





# ##########################################
# ##########################################
#FISH OUT V-GENE NT SEQUENCES FROM MERGED FILE
# targetClones <- read.csv("/Users/mfilipav/seqdata/imgt_align/cdr3_aa_input/march 14/mix.txt",
#                          header = FALSE, sep = "")
#
targetClones <- as.character(clones_for_fasta)
targetClones
targetClones[1]

sequences <- vector('list', length(targetClones))
length(targetClones)

for(i in 1:length(targetClones)){

    hMask <- which(imgt_5aa_prod_filtered$CDR3.IMGT == targetClones[i])
    sequences[[i]] <- unlist(strsplit(names(sort(table(imgt_5aa_prod_filtered$Sequence[hMask]), decreasing = TRUE)[1]), "*"))
}
getwd()
write.fasta(sequences, targetClones, 'bcma_vgenes.fasta')
# 

