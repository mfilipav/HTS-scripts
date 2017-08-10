# Last update: 2017-08-10
# 2nd BCR sequencing analysis pipeline using IMGT alignment output
# data obtained from SF's stand-alone maf aligner
# samples prepared using non-degenerate primer sets

# Analysis done using:
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
setwd("~/seqdata/2017-08-05/imgt_alignment/")

file_list = list.files(pattern = "*.txt")
#file_list <- file_list[1] # subset desired files
file_list

count.fields("data.txt", sep = "\t")

colclasses <- rep("NULL", 31)  # 31 cols in IMGT's 1_Summary file
colclasses[c(1, 3, 4, 5, 6)] <- 
    rep("character", length(c(1, 3, 4, 5, 6)), colclasses) # 5 and 31 too



# Pull data into a single list of DFs
data_list_dirty = lapply(file_list, read.table, fill = TRUE, header = TRUE,
                         sep = "\t", na.strings = " ",
                         colClasses = colclasses)
#View(data_list_dirty[[1]])


# Rename samples
# sample_names <- c("VH1a", 
#                   "VH1b", 
#                   "VH1c",
#                   "VH1d",
#                   "VH2a",
#                   "VH2b",
#                   "VH3a",
#                   "VH3b",
#                   "VH3c",
#                   "VH4",
#                   #"VH5a",
#                   "VH6a")
# names(data_list_dirty) <- sample_names



########### DATA CLEAN-UP ########
# Collect stats about raw dataset
data_stats <- list()
rm(data_stats)

for (i in length(data_list_dirty)) {
    data_stats[[1]] <- names(data_list_dirty[1]) # names
    data_stats[[2]] <- length(data_list_dirty[[1]]$Functionality)    
    data_stats[[3]] <- length(which(data_list_dirty[[1]]$Functionality == "no result"))    # not aligned
    # data_stats[i][3] <- length(which(data_list_dirty[[i]]$Functionality == "unproductive")) # aligned
}

View(data_list_dirty[[1]])


######## Filter for unproductive V-genes ##########
# Here we have only "unproductive" and "no result" labels
# None of the V-genes are productive, because of CDR3 missing
data_list <- lapply(data_list_dirty, function(df) {
    subset(df, df$V.DOMAIN.Functionality == "unknown (see comment)")
})
# Alt solution using which:
# data_list[[1]] <- data_list_dirty[[1]][which(data_list_dirty[[1]]$Functionality == "productive"), ]
rm(data_list_dirty)


####### Fix V-gene names ########
# Check all possible factor levels in categorical variable vector
V_gene_factors <- factor(data_list[[1]]$V.GENE.and.allele)
V_gene_factors[[1]]

# cleanup V-D-J gene names

data_list <- lapply(data_list, function(df) {
    df$V.GENE.and.allele <- gsub("Homsap ", "", df$V.GENE.and.allele);
    # matches "homosap" and removes it
    return(df)
})

data_list <- lapply(data_list, function(df) {
    df$V.GENE.and.allele <- gsub("S", "-", df$V.GENE.and.allele);
    return(df)
})

data_list <- lapply(data_list, function(df) {
    df$V.GENE.and.allele <- gsub("\\*.*", "", df$V.GENE.and.allele);
    return(df)
})

data_list <- lapply(data_list, function(df) {
    df$V.GENE.and.allele.short <- gsub("\\-.*", "", df$V.GENE.and.allele);
    # matches everything after "-" and removes it
    return(df)
})


data_list[[1]]$V.GENE.and.allele


##### Stats about sequence quality #######
# print df lengths BEFORE singleton removal
lapply(data_list_dirty, function(df) {
    cat(paste(length(df$Functionality)), "\n")
})

# print df lengths after singleton removal
lapply(data_list, function(df) {
    cat(paste(length(df$Functionality)), "\n")
})





########## PLOTS ##########

########## VDJ gene count and allele frequency plots ##########

##### Allele frequency plot #####
# Frequency rank V-gene occurence, and store in list of DFs
vgene_prop_table_sorted <- lapply(data_list, function(df) {
    data.frame(sort(prop.table(table(df$V.GENE.and.allele))*100, decreasing = TRUE))
})

# Filter out sequences with Freq > 1%
vgene_prop_table_sorted_1 <- lapply(vgene_prop_table_sorted, function(df) {
    df[which(df$Freq > 0.01),]
})



vgene_short_prop_table_sorted <- lapply(data_list, function(df) {
    data.frame(sort(prop.table(table(df$V.GENE.and.allele.short))*100, decreasing = TRUE))
})


# Print individual frequency graphs for each sample
# freqencies for each df add up to 1
# https://stackoverflow.com/questions/17368223/ggplot2-multi-group-histogram-with-in-group-proportions-rather-than-frequency
for (i in 1:length(vgene_prop_table_sorted_1)) {
    
    #bla <- nrow(vgene_prop_table_sorted[[i]])
    #i = 3
    vgene_plot_freq <- ggplot(bind_rows(vgene_prop_table_sorted_1[[i]], .id="df"),
                              aes(x = Var1, y = Freq, fill = df)) +

        geom_bar(stat="identity", position='dodge') +
        theme(axis.text.x=element_text(angle = +45, hjust = 1)) +
        labs(x="V-gene", y="Frequency (%)")

    # vgene_plot_freq

    ggsave(paste("0", i, "_vgene1percent.png", sep = ""), vgene_plot_freq)

}









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
# data_list23 <- data_list[1]
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

intersection_matrix <- matrix(, length(data_list23), length(data_list23))
for (i in 1:length(data_list23)) {
    
    for (j in 1:length(data_list23)) {
        
        len_intersect <- length(intersect(unique(data_list23[[i]]$aaSeqCDR3), 
                                          unique(data_list23[[j]]$aaSeqCDR3)))
        
        len_min <- min(length(unique(data_list23[[i]]$aaSeqCDR3)), 
                       length(unique(data_list23[[j]]$aaSeqCDR3)))
        
        intersection_matrix[i, j] <- len_intersect / len_min * 100
        
    }
}

intersection_matrix_pretty <- round(intersection_matrix, 1)
rownames(intersection_matrix_pretty) <- c("n1a", "n1b", "m2a", "m2b", 
                                          "n2", "m3a", "m3b", "n3")

colnames(intersection_matrix_pretty) <- c("n1a", "n1b", "m2a", "m2b", 
                                          "n2", "m3a", "m3b", "n3")

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

barplo
mean(intersection_oneVsAll)

?barplot

View(data_list)






# create graphics object with 3 libraries
venn_cdr3_overlap <- venn.diagram(list("1 Library" = A, "2 Library" = B),
                                  filename = NULL, print.mode = "raw",
                                  #set print.mode = "perc", for percentage display
                                  
                                  cex = 1, fontfamily = "sans",
                                  cat.cex = 1, cat.dist = c(0.1,0.1),
                                  cat.fontface = "bold",
                                  cat.fontfamily = "sans")
grid.draw(venn_cdr3_overlap)







View(data_list23[1])

