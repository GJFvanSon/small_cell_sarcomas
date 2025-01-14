############################################################################################
#                  Creating copynumber plots over multiple samples                         #
#                                                                                          #
#                                   Gijs van Son                                           #          
#                                  January 8 2025                                          #
#                         g.j.f.vanson@prinsesmaximacentrum.nl                             #                            
#                                                                                          #
############################################################################################
########### setting the environment #############
library(pheatmap)
library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(gtools)
library(grid)

setwd(dir = "/Users/gijsvanson/OneDrive - Prinses Maxima Centrum/Femke/CIC_DUX4/")

##### read some ratio files #########

filenames <- list.files('DNA_analysis/CNV_files/')

ratio_files_list <- lapply(filenames[c(2,4,6,10, 17, 31, 33,19,45,22,44,29 )], function(x){
  fread(paste0("DNA_analysis/CNV_files/",x))
}
)
# make the names compatible with the rest of the paper
names(ratio_files_list) <-c("ES-016_TO", "ES-032_TO", "ES-055_TO", "ES-080_TO",
                            "ES-016_TT", "ES-055_TT", "ES-080_TT", "ES-090_TT",
                            "ES-090_TO", "ES-032_TT", "ES-045_TO", "ES-045_TT")

# put all the ratio files in one big table
ratio_table <- melt(ratio_files_list, id.vars = c("Chromosome", "Start", "Ratio", "MedianRatio", "CopyNumber"))

#set some variables (binsize changes for WES or WGS data)
df <- data.frame()
ploidy <- 2
binsize = 200000
maxLevelToPlot <- 4
chromosomes <- mixedsort(as.character(unique((ratio_table$Chromosome))))

# make new tables for all the samples to check gains and losses
for (sample in names(ratio_files_list)) {
  message(sample)
  ratio<- ratio_files_list[[sample]]
  ratio <- subset(ratio, Ratio!=-1)
  ratio <- ratio[!is.na(ratio$Ratio),]
  for (i in c(1:length(chromosomes))) {
    chrom <- chromosomes[i]
    message(chrom)
    # Select chromosome
    tmp <- subset( ratio, Chromosome==chrom)
    
    # Binning of positions
    grouping <- cut(tmp$Start, c(seq(0,max(ratio$Start[which(ratio$Chromosome == chrom)]),binsize), max(ratio$Start[which(ratio$Chromosome == chrom)])), 
                    labels=seq(binsize/2, (max(ratio$Start[which(ratio$Chromosome == chrom)])+(binsize/2)), binsize))
    # Calculate median ratio
    medians <- data.frame(tapply(tmp$Ratio, grouping, median))
    
    colnames(medians) <- "Ratio"
    # Alternate chromosome banding
    if ((i%%2) == 0) {
      medians$Cols <- 'neutral1'
    } else {
      medians$Cols <- 'neutral2'	
    }
    
    # Create clean data farme
    medians$Chromosome <- chrom
    medians$GenomicPosition <- as.numeric(as.character(rownames(medians)))
    medians$CopyNumber <- medians$Ratio*ploidy
    
    # Add event coloring
    medians$Cols[medians$CopyNumber>=ploidy+0.8] <- "gain"
    medians$Cols[medians$CopyNumber<=ploidy-0.8] <- "loss"
    
    # Merge data
    medians$Sample <- sample
    df <- rbind(df, medians)
    
  }
  df$CopyNumber[df$CopyNumber > maxLevelToPlot] <- maxLevelToPlot
  df$Chromosome <- factor(df$Chromosome, levels=chromosomes)
}
# put the data in a new frame
ratio_frame <- data.table()
ratio_frame$Chromosome <- df$Chromosome[which(df$Sample == names(ratio_files_list)[1])]
ratio_frame$GenomicPosition <-  df$GenomicPosition[which(df$Sample == names(ratio_files_list)[1])]
# fill the frame for each sample
for(sample in names(ratio_files_list)){
  message(sample)
  ratio_frame$temp_ratio <- NA
  ratio_frame$temp_ratio[which(paste0(ratio_frame$Chromosome,
                                      "_", ratio_frame$GenomicPosition) %in% 
                                 paste0(df$Chromosome[which(df$Sample == sample)],
                                        "_", df$GenomicPosition[which(df$Sample == sample)]))] <- 
    df$Ratio[intersect(which(df$Sample == sample),
                       which(paste0(df$Chromosome, "_", df$GenomicPosition) %in% 
                               paste0(ratio_frame$Chromosome, "_", ratio_frame$GenomicPosition)))]
  ratio_frame$temp_copy <- NA
  ratio_frame$temp_copy[which(paste0(ratio_frame$Chromosome, 
                                     "_", ratio_frame$GenomicPosition) %in% 
                                paste0(df$Chromosome[which(df$Sample == sample)],
                                       "_", df$GenomicPosition[which(df$Sample == sample)]))] <-
    df$CopyNumber[intersect(which(df$Sample == sample),
                            which(paste0(df$Chromosome, "_", df$GenomicPosition) %in% paste0(ratio_frame$Chromosome, "_", ratio_frame$GenomicPosition)))]
  colnames(ratio_frame)[which(colnames(ratio_frame) == "temp_ratio")] <- paste0(sample, "_Ratio")
  colnames(ratio_frame)[which(colnames(ratio_frame) == "temp_copy")] <- paste0(sample, "_Copynumber")
}

ratio_frame <- as.data.frame(ratio_frame)
# find all NA values in the dataframe and remove those lines
cleaning_list <- lapply(colnames(ratio_frame)[-c(1,2)], function(x){
  which(is.na(ratio_frame[,x]))
})
cleaning_list <- unique(unlist(cleaning_list))
ratio_frame_clean <- ratio_frame[-cleaning_list,]
# make a heatmap for all the samples over the chromosome
heatmap_table_copynumber <- t(ratio_frame_clean[,grep("Copynumber", colnames(ratio_frame_clean))])
annotation_df <- data.frame("Chromosome" = ratio_frame_clean[c(which(ratio_frame_clean$Chromosome == 1),
                                                               which(ratio_frame_clean$Chromosome == 2),
                                                               which(ratio_frame_clean$Chromosome == 3),
                                                               which(ratio_frame_clean$Chromosome == 4),
                                                               which(ratio_frame_clean$Chromosome == 5),
                                                               which(ratio_frame_clean$Chromosome == 6),
                                                               which(ratio_frame_clean$Chromosome == 7),
                                                               which(ratio_frame_clean$Chromosome == 8),
                                                               which(ratio_frame_clean$Chromosome == 9),
                                                               which(ratio_frame_clean$Chromosome == 10),
                                                               which(ratio_frame_clean$Chromosome == 11),
                                                               which(ratio_frame_clean$Chromosome == 12),
                                                               which(ratio_frame_clean$Chromosome == 13),
                                                               which(ratio_frame_clean$Chromosome == 14),
                                                               which(ratio_frame_clean$Chromosome == 15),
                                                               which(ratio_frame_clean$Chromosome == 16),
                                                               which(ratio_frame_clean$Chromosome == 17),
                                                               which(ratio_frame_clean$Chromosome == 18),
                                                               which(ratio_frame_clean$Chromosome == 19),
                                                               which(ratio_frame_clean$Chromosome == 20),
                                                               which(ratio_frame_clean$Chromosome == 21),
                                                               which(ratio_frame_clean$Chromosome == 22)),1])
rownames(annotation_df) <- colnames(heatmap_table_copynumber[,c(which(ratio_frame_clean$Chromosome == 1),
                                                                                                                   which(ratio_frame_clean$Chromosome == 2),
                                                                                                                   which(ratio_frame_clean$Chromosome == 3),
                                                                                                                   which(ratio_frame_clean$Chromosome == 4),
                                                                                                                   which(ratio_frame_clean$Chromosome == 5),
                                                                                                                   which(ratio_frame_clean$Chromosome == 6),
                                                                                                                   which(ratio_frame_clean$Chromosome == 7),
                                                                                                                   which(ratio_frame_clean$Chromosome == 8),
                                                                                                                   which(ratio_frame_clean$Chromosome == 9),
                                                                                                                   which(ratio_frame_clean$Chromosome == 10),
                                                                                                                   which(ratio_frame_clean$Chromosome == 11),
                                                                                                                   which(ratio_frame_clean$Chromosome == 12),
                                                                                                                   which(ratio_frame_clean$Chromosome == 13),
                                                                                                                   which(ratio_frame_clean$Chromosome == 14),
                                                                                                                   which(ratio_frame_clean$Chromosome == 15),
                                                                                                                   which(ratio_frame_clean$Chromosome == 16),
                                                                                                                   which(ratio_frame_clean$Chromosome == 17),
                                                                                                                   which(ratio_frame_clean$Chromosome == 18),
                                                                                                                   which(ratio_frame_clean$Chromosome == 19),
                                                                                                                   which(ratio_frame_clean$Chromosome == 20),
                                                                                                                   which(ratio_frame_clean$Chromosome == 21),
                                                                                                                   which(ratio_frame_clean$Chromosome == 22))])
heatmap_table_copynumber <- t(ratio_frame_clean[,grep("Copynumber", colnames(ratio_frame_clean))])
# set the names for the figure
names4fig <- rownames(heatmap_table_copynumber)
table4fig <- heatmap_table_copynumber[which(rownames(heatmap_table_copynumber) %in% names4fig),c(which(ratio_frame_clean$Chromosome == 1),
                                                                                                 which(ratio_frame_clean$Chromosome == 2),
                                                                                                 which(ratio_frame_clean$Chromosome == 3),
                                                                                                 which(ratio_frame_clean$Chromosome == 4),
                                                                                                 which(ratio_frame_clean$Chromosome == 5),
                                                                                                 which(ratio_frame_clean$Chromosome == 6),
                                                                                                 which(ratio_frame_clean$Chromosome == 7),
                                                                                                 which(ratio_frame_clean$Chromosome == 8),
                                                                                                 which(ratio_frame_clean$Chromosome == 9),
                                                                                                 which(ratio_frame_clean$Chromosome == 10),
                                                                                                 which(ratio_frame_clean$Chromosome == 11),
                                                                                                 which(ratio_frame_clean$Chromosome == 12),
                                                                                                 which(ratio_frame_clean$Chromosome == 13),
                                                                                                 which(ratio_frame_clean$Chromosome == 14),
                                                                                                 which(ratio_frame_clean$Chromosome == 15),
                                                                                                 which(ratio_frame_clean$Chromosome == 16),
                                                                                                 which(ratio_frame_clean$Chromosome == 17),
                                                                                                 which(ratio_frame_clean$Chromosome == 18),
                                                                                                 which(ratio_frame_clean$Chromosome == 19),
                                                                                                 which(ratio_frame_clean$Chromosome == 20),
                                                                                                 which(ratio_frame_clean$Chromosome == 21),
                                                                                                 which(ratio_frame_clean$Chromosome == 22))]
table4fig <- table4fig[order(rownames(table4fig)),]
# plot the figure
pheatmap(table4fig,
         cluster_cols = F, 
         show_colnames = F, 
         annotation_col = data.frame("Chromosome" = as.character(ratio_frame_clean$Chromosome[c(which(ratio_frame_clean$Chromosome == 1),
                                                                                                which(ratio_frame_clean$Chromosome == 2),
                                                                                                which(ratio_frame_clean$Chromosome == 3),
                                                                                                which(ratio_frame_clean$Chromosome == 4),
                                                                                                which(ratio_frame_clean$Chromosome == 5),
                                                                                                which(ratio_frame_clean$Chromosome == 6),
                                                                                                which(ratio_frame_clean$Chromosome == 7),
                                                                                                which(ratio_frame_clean$Chromosome == 8),
                                                                                                which(ratio_frame_clean$Chromosome == 9),
                                                                                                which(ratio_frame_clean$Chromosome == 10),
                                                                                                which(ratio_frame_clean$Chromosome == 11),
                                                                                                which(ratio_frame_clean$Chromosome == 12),
                                                                                                which(ratio_frame_clean$Chromosome == 13),
                                                                                                which(ratio_frame_clean$Chromosome == 14),
                                                                                                which(ratio_frame_clean$Chromosome == 15),
                                                                                                which(ratio_frame_clean$Chromosome == 16),
                                                                                                which(ratio_frame_clean$Chromosome == 17),
                                                                                                which(ratio_frame_clean$Chromosome == 18),
                                                                                                which(ratio_frame_clean$Chromosome == 19),
                                                                                                which(ratio_frame_clean$Chromosome == 20),
                                                                                                which(ratio_frame_clean$Chromosome == 21),
                                                                                                which(ratio_frame_clean$Chromosome == 22))]), 
                                     row.names = rownames(ratio_frame_clean)[c(which(ratio_frame_clean$Chromosome == 1),
                                                                               which(ratio_frame_clean$Chromosome == 2),
                                                                               which(ratio_frame_clean$Chromosome == 3),
                                                                               which(ratio_frame_clean$Chromosome == 4),
                                                                               which(ratio_frame_clean$Chromosome == 5),
                                                                               which(ratio_frame_clean$Chromosome == 6),
                                                                               which(ratio_frame_clean$Chromosome == 7),
                                                                               which(ratio_frame_clean$Chromosome == 8),
                                                                               which(ratio_frame_clean$Chromosome == 9),
                                                                               which(ratio_frame_clean$Chromosome == 10),
                                                                               which(ratio_frame_clean$Chromosome == 11),
                                                                               which(ratio_frame_clean$Chromosome == 12),
                                                                               which(ratio_frame_clean$Chromosome == 13),
                                                                               which(ratio_frame_clean$Chromosome == 14),
                                                                               which(ratio_frame_clean$Chromosome == 15),
                                                                               which(ratio_frame_clean$Chromosome == 16),
                                                                               which(ratio_frame_clean$Chromosome == 17),
                                                                               which(ratio_frame_clean$Chromosome == 18),
                                                                               which(ratio_frame_clean$Chromosome == 19),
                                                                               which(ratio_frame_clean$Chromosome == 20),
                                                                               which(ratio_frame_clean$Chromosome == 21),
                                                                               which(ratio_frame_clean$Chromosome == 22))]),
         annotation_legend = F,
         annotation_names_col = T, 
         annotation_names_row = T,
         main = "Copynumber of organoid samples",
         gaps_col =  unlist(lapply(c(1:22), function(x){grep(x, ratio_frame_clean$Chromosome)[1]-1})),
         cluster_rows = F,
         color = c("darkblue", "white", "red"),
         annotation_colors = list(Chromosome = c('1' = "black",'2'="gray",'3'= "black",'4'="gray",'5'="black",'6'="gray", '7'="black", '8'="gray",
                                                 '9'="black",'10'="gray",'11'= "black",'12'="gray",'13'="black", '14'="gray", '15'="black", '16'="gray",
                                                 '17'="black",'18'= "gray",'19'= "black",'20'= "gray", '21'="black", '22'="gray"))
)
# add the chromosome numbers (when possible, also possible to do this later)
grid.text(as.factor(unique(annotation_df$Chromosome)), x=c(0.05,0.125,0.185,0.238,0.28,0.33,0.38,0.42,0.46,0.50,
                                                           0.545,0.59,0.623,0.65,0.68,0.713,0.745,0.772, 0.80,0.827,0.848,0.863),y=rep(0.89,22), gp=gpar(fontsize=12, col = "white", fontface = "bold"))
dev.off()



#making a copynumber plot overlaying the different samples from sampleset es16.
df_tissue <- data.frame("Chromosome" = ratio_frame_clean$Chromosome,
                        "GenomicPosition" = ratio_frame_clean$GenomicPosition,
                        "CopyNumber" = ratio_frame_clean$`ES-016_TT_Copynumber`,
                        "Sample" = "Tissue"
)               
df_organoid <- data.frame("Chromosome" = ratio_frame_clean$Chromosome,
                          "GenomicPosition" = ratio_frame_clean$GenomicPosition,
                          "CopyNumber" = ratio_frame_clean$`ES-016_TO_Copynumber`,
                          "Sample" = "Organoid"
)
df <- rbind(df_tissue, df_organoid)
p <- ggplot(df, aes(GenomicPosition, CopyNumber, group=Chromosome, alpha=0.5)) +
  geom_point(aes(colour=Sample), shape=20, size=0.2, alpha=0.2) + 
  ylim(0,maxLevelToPlot) + 
  facet_wrap(~Chromosome, nrow=1, scales="free_x") + 
  scale_colour_manual(name="Sample", values=c("steelblue", "firebrick")) 

# Make it look clean
p <- p + theme(
  axis.line.x=element_blank(),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  panel.grid.minor=element_blank(),
  panel.grid.major.x=element_blank(),
  panel.grid.major.y=element_line(colour="grey80", linetype=2, size=0.2),
  legend.position="right",
  panel.background=element_blank(),
  panel.border=element_blank(),
  plot.background=element_blank())
print(p)

#making a copynumber plot overlaying the different samples from sampleset es32.
df_tissue <- data.frame("Chromosome" = ratio_frame_clean$Chromosome,
                        "GenomicPosition" = ratio_frame_clean$GenomicPosition,
                        "CopyNumber" = ratio_frame_clean$`ES-032_TT_Copynumber`,
                        "Sample" = "Tissue"
)               
df_organoid <- data.frame("Chromosome" = ratio_frame_clean$Chromosome,
                          "GenomicPosition" = ratio_frame_clean$GenomicPosition,
                          "CopyNumber" = ratio_frame_clean$`ES-032_TO_Copynumber`,
                          "Sample" = "Organoid"
)
df <- rbind(df_tissue, df_organoid)
p <- ggplot(df, aes(GenomicPosition, CopyNumber, group=Chromosome, alpha=0.5)) +
  geom_point(aes(colour=Sample), shape=20, size=0.2, alpha=0.2) + 
  ylim(0,maxLevelToPlot) + 
  facet_wrap(~Chromosome, nrow=1, scales="free_x") + 
  scale_colour_manual(name="Sample", values=c("steelblue", "firebrick")) 

# Make it look clean
p <- p + theme(
  axis.line.x=element_blank(),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  panel.grid.minor=element_blank(),
  panel.grid.major.x=element_blank(),
  panel.grid.major.y=element_line(colour="grey80", linetype=2, size=0.2),
  legend.position="right",
  panel.background=element_blank(),
  panel.border=element_blank(),
  plot.background=element_blank())
print(p)


