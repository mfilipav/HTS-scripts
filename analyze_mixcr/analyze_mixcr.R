# Last update: 2017-08-12
# 2nd BCR sequencing analysis pipeline using MiXCR stand alone alignment
# for file in *.fasta; do mixcr align -OallowPartialAlignments=true 
# --chains IGH --species hsa --report alignment_"${file}".log --save-description 
# "${file}" alignment_"${file}".vdjca; done
# 
# for file in *.vdjca; do mixcr exportAlignments --no-spaces -vGene -vHitScore 
# "${file}" alignment_"${file}".txt; done

# MiXCR v2.1.2 (built Thu Apr 20 11:22:25 CEST 2017; rev=6abb137; branch=hotfix/v2.1.2)
rm(list=ls())
getwd()
library(ggplot2)
library(gplots) # for heatmap.2
library(seqinr)
library(tidyverse)
library(VennDiagram)
library(reshape2)
library(data.table)
# library(NMF)              # heatmap, gaujoux and seoighe, 2010
# library(ComplexHeatmap)   # heatmap, gu, 2015


###### Import data #####
setwd("~/seqdata/2017-08-21/mixcr/clones/")

file_list = list.files(pattern = "*txt")
#file_list <- file_list[1] # subset desired files
file_list

count.fields("data.txt", sep = "\t")

colclasses <- rep("NULL", 35)  # 31 cols in IMGT's 1_Summary file
colclasses[c(2, 3, 6, 33)] <- 
    rep("character", length(c(2, 3, 6, 33)), colclasses) # 5 and 31 too

# Pull data into a single list of DFs
data_list_dirty = lapply(file_list, read.table, fill = TRUE, header = TRUE,
                         sep = "\t", na.strings = " ",
                         colClasses = colclasses)
View(data_list_dirty[[1]])
#write.csv(data_list_df, "./clonotpyes_12pcr_products.csv", sep="\t")


# Rename samples
sample_names <- c("VH1a",
                  "VH1b",
                  "VH1c",
                  "VH1d",
                  "VH2a",
                  "VH2b",
                  "VH3a",
                  "VH3b",
                  "VH3c",
                  "VH4",
                  "VH5a",
                  "VH6a")
names(data_list_dirty) <- sample_names



########### DATA CLEAN-UP ########

######## Filter for unproductive V-genes ##########
# remove stop codons
parasites <- c("\\*", "\\_")
data_list <- lapply(data_list_dirty, function(df)
    df[!grepl(paste(parasites, collapse = "|"), df$aaSeqCDR3), ]
)

# cleanup V-D-J gene names
data_list <- lapply(data_list, function(df) {
    df$allVHitsWithScore <- gsub("\\*.*", "", df$allVHitsWithScore);
    #df$allDHitsWithScore <- gsub("\\*.*", "", df$allDHitsWithScore);
    #df$allJHitsWithScore <- gsub("\\*.*", "", df$allJHitsWithScore);
    return(df)
})

# add CDR3 nt and aa lengths as new df columns
data_list <- lapply(data_list, function(df) {
    #df$nSeqCDR3Length <- nchar(df$nSeqCDR3); 
    df$aaSeqCDR3Length <- nchar(df$aaSeqCDR3);
    return(df)
})


######## Filter out empty rows ##########
data_list <- lapply(data_list, function(df) {
    subset(df, df$allVHitsWithScore != "")
})

########### Adding ID after cleanup #############
data_list_with_singletons <- lapply(data_list, function(df) {
    df$id <- seq.int(nrow(df));
    return(df)
})

########### Singleton Removal #############
data_list <- lapply(data_list_with_singletons, function(df) {
    subset(df, df$cloneCount > 1)
})


########### Update/Normalize cloneFractions to % value#######
data_list <- lapply(data_list, function(df) {
    df$cloneFraction <- as.numeric(df$cloneFraction) / sum(as.numeric(df$cloneFraction)) * 100;
    return(df)
    # sum(as.numeric(test_data$cloneFractionNorm))  ## test sum for 100%
})




##### Stats about sequence quality #######
# print df lengths BEFORE singleton removal
lapply(data_list, function(df) {
    cat(paste(length(df$cloneCount)), "\n")
})

# print df lengths after singleton removal
lapply(data_list_wout_singletons, function(df) {
    cat(paste(length(df$cloneCount)), "\n")
})




########## PLOTS ##########
data_list_wsingletons <- data_list
data_list <- data_list_wout_singletons

########## VDJ gene count and allele frequency plots ##########

##### Allele frequency plot #####
# Frequency rank V-gene occurence, and store in list of DFs
vgene_prop_table_sorted <- lapply(data_list, function(df) {
    data.frame(sort(prop.table(table(df$allVHitsWithScore))*100, decreasing = TRUE))
})


vgene_freq_df <- bind_rows(vgene_prop_table_sorted, .id = "df" )
write.table(vgene_freq_df, "./vgene_frequencies_no_singletons_20170821_mf.txt", sep="\t")


# Filter out sequences with Freq > 1%
vgene_prop_table_sorted_1 <- lapply(vgene_prop_table_sorted, function(df) {
    df[which(df$Freq > 1),]
})



# Print individual frequency graphs for each sample
# freqencies for each df add up to 1
# https://stackoverflow.com/questions/17368223/ggplot2-multi-group-histogram-with-in-group-proportions-rather-than-frequency
for (i in 1:length(vgene_prop_table_sorted_1)) {
    
    #bla <- nrow(vgene_prop_table_sorted[[i]])
    #i = 3
    df = names(vgene_prop_table_sorted_1[i])
    vgene_plot_freq <- ggplot(bind_rows(vgene_prop_table_sorted_1[[i]], 
                                        .id=df),
                              aes(x = Var1, y = Freq, fill = df)) +

        geom_bar(stat="identity", position='dodge') +
        theme(axis.text.x=element_text(angle = +45, hjust = 1)) +
        labs(x="V-gene", y="Frequency (%)")

    # vgene_plot_freq

    ggsave(filename = paste("0", i, "_", df, "_nosingletons_vgene_1percent.pdf", sep = ""), 
           plot = vgene_plot_freq, width = 8.97, height = 8.97, device = "pdf")

}


################### 1 2 3  plot ####################
data_list_1 <- data_list[5:6]
data_list_2 <- data_list[c(1, 4, 8)]
data_list_3 <- data_list[c(2, 3, 7)]

for (i in 1:length(vgene_prop_table_sorted_1)) {
    
    #bla <- nrow(vgene_prop_table_sorted[[i]])
    #i = 3
    df = names(vgene_prop_table_sorted_1[i])
    vgene_plot_freq <- ggplot(bind_rows(vgene_prop_table_sorted_1[[i]], 
                                        .id=df),
                              aes(x = Var1, y = Freq, fill = df)) +
        
        geom_bar(stat="identity", position='dodge') +
        theme(axis.text.x=element_text(angle = +45, hjust = 1)) +
        labs(x="V-gene", y="Frequency (%)")
    
    # vgene_plot_freq
    
    ggsave(filename = paste("0", i, "_", df, "sample_1_vgene_1percent.pdf", sep = ""), 
           plot = vgene_plot_freq, device = "pdf")
}


## 2017-08-22
########## Species Accumulation Plots #######

test_data <- data_list[[1]]


# Add id indeces {1, 2, 3, ... , nrow(df)}
test_data$id <- seq.int(nrow(test_data))
cdr3_all <- test_data$aaSeqCDR3
cdr3_unique <- length(unique(cdr3_all))

## print unique clones for each data set:

lapply(data_list, function(df) {
    cat(paste(length(unique(df$aaSeqCDR3))), "\n")
})

sum(as.numeric(test_data$cloneFraction)) 

# normalize cloneFractions to account for cleanup
test_data$cloneFractionNorm <- as.numeric(test_data$cloneFraction) / sum(as.numeric(test_data$cloneFraction)) * 100
# sum(as.numeric(test_data$cloneFractionNorm))  ## test sum for 100%

for (i in 1:length(vgene_prop_table_sorted_1)) {
    
    #bla <- nrow(vgene_prop_table_sorted[[i]])
    #i = 3
    #df = names(vgene_prop_table_sorted_1[i])
    
    cdr3_frequency_plot <- ggplot(test_data, aes(x = id, y = cloneFractionNorm)) + 
                            geom_point(colour = "black") +
                            labs(x="Clone count", y="CDR3 Clone Frequency (%)")
    
    ggsave(filename = paste("0", i, "_", df, "cdr3_freq.pdf", sep = ""), 
           plot = cdr3_frequency_plot, device = "pdf")
}

rm(cdr3_frequency_plot)











########## Barplot number clones ##########

# Barplot with number of clones
length_clones <- as.list(length)
names(length_clones) <- sample_names
?as.list
barplot(length_clones, ylab = "Number of Clones")

barplot(1:12, xlab = "bla")



# # plots v-gene counts dodged for multiple samples
# ggplot(bind_rows(data_list, .id="df")) +
#     geom_bar(
#         mapping = aes(V.GENE.and.allele.short, fill = df),
#         position = "dodge") +
#     theme(axis.text.x=element_text(angle = +90, hjust = 0)) # vertical axis

# # for individual samples percent frequency is calculated among all samples 
# ggplot(bind_rows(data_list[[1]]), aes(x = factor(V.GENE.and.allele.short))) +  
#     geom_bar(aes(y = (..count..)/sum(..count..))) + 
#     scale_y_continuous(labels = scales::percent) +
#     theme(axis.text.x=element_text(angle = +45, hjust = 1)) +
#     labs(x="", y="Frequency (%)")
# 

# # percentage plots 
# data_list <- data_list[1]
# ggplot(bind_rows(data_list23, .id="df"), aes(x = factor(V.GENE.and.allele), fill = df)) +  
#     geom_bar(aes(y = (..count..)/sum(..count..)), position = "dodge") + 
#     scale_y_continuous(labels = scales::percent) +
#     theme(axis.text.x=element_text(angle = +45, hjust = 1)) +
#     labs(x="", y="Frequency (%)")
# 
# 
# # do percent frequencies for each samples separately!!
# data_list23 <- data_list[5:6]
# ggplot(bind_rows(data_list23, .id="df"), aes(x = factor(allVHitsWithScore), fill = df)) +  
#     geom_bar(aes(y = (..count..)/sum(..count..)), position = "dodge") + 
#     scale_y_continuous(labels = scales::percent) +
#     theme(axis.text.x=element_text(angle = +45, hjust = 1)) +
#     labs(x="", y="Frequency (%)")
# 

# # count barplots
# vgene_plot_count <- ggplot(bind_rows(data_list, .id="df"), aes(x=V.GENE.and.allele.short, fill = df)) + 
#     geom_histogram(binwidth=.5, alpha=.5, position="dodge", stat="count")
# 
# vgene_plot_freq <- ggplot(bind_rows(data_list, .id="df"), aes(x=V.GENE.and.allele.short)) + 
#     geom_bar(aes(y= (..count..)/sum(..count..), fill = df, group = df), position="dodge")
# 
# 
# theme_bw() + geom_bar(stat="identity", fill = df, width=0.5) +
#     labs(x="", y="Frequency (%)") +
#     theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1)) +
#     theme(axis.title.y = element_text(colour="black", size=15, vjust=1.5, hjust=0.4),
#           axis.text.y = element_text(angle=0, vjust=0.5, size=10))
# 
# ## 




########## Rank Plot #######
rank_plots <- as.list(0)
cdr3_freq <- data.frame(as.numeric(data_list[[1]]$cloneFraction)*10000)
cdr3_freq$index <- 1:nrow(cdr3_freq)
colnames(cdr3_freq)[1] <- "V1"
cdr3_freq <- rbind(cdr3_freq, data.frame(V1=cumsum(cdr3_freq$V1),
                                         index=1:nrow(cdr3_freq)))

cdr3_freq$facet <- factor(rep(factor(c("Frequency","Cumulative frequency")),
                              each=length(data_list[[1]]$cloneFraction)),
                          levels = c("Frequency", "Cumulative frequency"))

rank_plots[[1]] <- ggplot(cdr3_freq, aes(x=index, y= V1))
rank_plots[[1]] <- rank_plots[[1]] + geom_point(size=1, alpha = 0.8) +
    facet_grid(facet ~., scales = "free") +
    labs(x= "CDR3 index", y="CDR3 frequency (%)") +
    theme_bw() +
    theme(plot.title = element_text(face="bold", size=rel(1), hjust=0),
          plot.margin = unit(c(2, 2, 2, 2),"points"),
          strip.text = element_text(size = rel(0.75)),
          legend.position = "none") +
    theme(strip.background = element_rect(fill = scales:::alpha("blue", 0.3)))
rank_plots




#### OVERLAP AND INTERSECTION %
# must figure out why sequences that are duplicated don't appear so
# if "duplicated" are not removed, or alternatively, if "unique" are not kept,
# then I won't have 100% identity across the columns!!!
data_list23 <- data_list[1:8]

intersection_matrix <- matrix(, length(data_list), length(data_list))
for (i in 1:length(data_list)) {
    
    for (j in 1:length(data_list)) {
        
        len_intersect <- length(intersect(unique(data_list[[i]]$aaSeqCDR3), 
                                          unique(data_list[[j]]$aaSeqCDR3)))
        
        len_min <- min(length(unique(data_list[[i]]$aaSeqCDR3)), 
                       length(unique(data_list[[j]]$aaSeqCDR3)))
        
        intersection_matrix[i, j] <- len_intersect / len_min * 100
        
    }
}

intersection_matrix_pretty <- round(intersection_matrix, 1)
rownames(intersection_matrix_pretty) <- sample_names

colnames(intersection_matrix_pretty) <- sample_names

intersection_matrix_pretty
col_scheme <- colorRampPalette(c("yellow", "red"))(n = 1000)
heatmap.2(intersection_matrix_pretty, Colv = FALSE, Rowv = FALSE, col = col_scheme, tracecol="black")



# Melt list into df, number of unique clones
melt_data_list <- melt(data_list)
length((melt_data_list$aaSeqCDR3))
length(unique(melt_data_list$aaSeqCDR3))



# Show top 100 clones for each 12 samples, how similar are they to each other?
top10cdr3 <- lapply(data_list, function(df) {
    df$aaSeqCDR3[1:1500]
})


melt_top10cdr3 <- melt(top10cdr3)

unique(melt_top10cdr3$value)


data_list <- data_list[1:8]

# Overlap between whole sample set and individual samples
intersection_oneVsAll <- as.numeric()
for (i in 1:length(data_list)) {
    
    compare_against <- melt_data_list[which(melt_data_list$L1 != i), ]
    len_intersect <- length(intersect(unique(data_list[[i]]$aaSeqCDR3), 
                                      unique(compare_against$aaSeqCDR3)))
    
    len_min <- length(unique(data_list[[i]]$aaSeqCDR3))
    intersection_oneVsAll[i] <- len_intersect / len_min * 100
}
xlabs <- c("2a", "2b", "2c", "3a", "3b", "3c")
barplot(intersection_oneVsAll,
        ylab = "CDR3 aa overlap (%) between individual and all samples",
)



