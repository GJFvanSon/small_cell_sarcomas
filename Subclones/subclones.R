############################################################################################
#                     DNA subclonal analysis on Ewing(-like) sarcomas                      #
#                                                                                          #
#                                   Gijs van Son                                           #          
#                                   June 11 2024                                           #
#                         g.j.f.vanson@prinsesmaximacentrum.nl                             #                            
#                                                                                          #
############################################################################################
########### setting the environment #############
library(ggplot2)
library(Matrix)
library(ggalluvial)
library(reshape2)
library(Seurat)
setwd("/Users/gijsvanson/OneDrive - Prinses Maxima Centrum/Femke/CIC_DUX4/DNA_analysis/Subclones/")

# read the output form pyclone-vi (a program to identify subclones)
pyclone_out <- read.table("pyclone_output.tsv", header = T)

# lets have a look at whats in there
picture <- ggplot(pyclone_out) + geom_point(aes(size = cellular_prevalence, x = sample_id, y = as.vector(cluster_id), color=as.character(cluster_id)))
############# data wrangling to get the right input format ###########
# construct a frequency table for later figures
fraq_table <- table(pyclone_out[,c("sample_id","cluster_id")])
fraq_table["ES16O", "0"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16O" & pyclone_out$cluster_id == 0), 4])
fraq_table["ES16O", "1"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16O" & pyclone_out$cluster_id == 1), 4])
fraq_table["ES16O", "2"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16O" & pyclone_out$cluster_id == 2), 4])
fraq_table["ES16O", "3"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16O" & pyclone_out$cluster_id == 3), 4])
fraq_table["ES16O", "4"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16O" & pyclone_out$cluster_id == 4), 4])
fraq_table["ES16O", "5"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16O" & pyclone_out$cluster_id == 5), 4])
fraq_table["ES16O", "6"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16O" & pyclone_out$cluster_id == 6), 4])
fraq_table["ES16O", "7"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16O" & pyclone_out$cluster_id == 7), 4])
fraq_table["ES16O", "8"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16O" & pyclone_out$cluster_id == 8), 4])
fraq_table["ES16O", "9"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16O" & pyclone_out$cluster_id == 9), 4])
fraq_table["ES16O", "10"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16O" & pyclone_out$cluster_id == 10), 4])
fraq_table["ES16O", "11"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16O" & pyclone_out$cluster_id == 11), 4])
fraq_table["ES16O", "12"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16O" & pyclone_out$cluster_id == 12), 4])
fraq_table["ES16O", "13"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16O" & pyclone_out$cluster_id == 13), 4])
fraq_table["ES16O", "14"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16O" & pyclone_out$cluster_id == 14), 4])
fraq_table["ES16O", "15"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16O" & pyclone_out$cluster_id == 15), 4])
fraq_table["ES16T", "0"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16T" & pyclone_out$cluster_id == 0), 4])
fraq_table["ES16T", "1"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16T" & pyclone_out$cluster_id == 1), 4])
fraq_table["ES16T", "2"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16T" & pyclone_out$cluster_id == 2), 4])
fraq_table["ES16T", "3"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16T" & pyclone_out$cluster_id == 3), 4])
fraq_table["ES16T", "4"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16T" & pyclone_out$cluster_id == 4), 4])
fraq_table["ES16T", "5"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16T" & pyclone_out$cluster_id == 5), 4])
fraq_table["ES16T", "6"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16T" & pyclone_out$cluster_id == 6), 4])
fraq_table["ES16T", "7"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16T" & pyclone_out$cluster_id == 7), 4])
fraq_table["ES16T", "8"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16T" & pyclone_out$cluster_id == 8), 4])
fraq_table["ES16T", "9"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16T" & pyclone_out$cluster_id == 9), 4])
fraq_table["ES16T", "10"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16T" & pyclone_out$cluster_id == 10), 4])
fraq_table["ES16T", "11"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16T" & pyclone_out$cluster_id == 11), 4])
fraq_table["ES16T", "12"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16T" & pyclone_out$cluster_id == 12), 4])
fraq_table["ES16T", "13"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16T" & pyclone_out$cluster_id == 13), 4])
fraq_table["ES16T", "14"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16T" & pyclone_out$cluster_id == 14), 4])
fraq_table["ES16T", "15"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES16T" & pyclone_out$cluster_id == 15), 4])
fraq_table["ES55O", "0"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55O" & pyclone_out$cluster_id == 0), 4])
fraq_table["ES55O", "1"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55O" & pyclone_out$cluster_id == 1), 4])
fraq_table["ES55O", "2"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55O" & pyclone_out$cluster_id == 2), 4])
fraq_table["ES55O", "3"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55O" & pyclone_out$cluster_id == 3), 4])
fraq_table["ES55O", "4"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55O" & pyclone_out$cluster_id == 4), 4])
fraq_table["ES55O", "5"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55O" & pyclone_out$cluster_id == 5), 4])
fraq_table["ES55O", "6"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55O" & pyclone_out$cluster_id == 6), 4])
fraq_table["ES55O", "7"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55O" & pyclone_out$cluster_id == 7), 4])
fraq_table["ES55O", "8"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55O" & pyclone_out$cluster_id == 8), 4])
fraq_table["ES55O", "9"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55O" & pyclone_out$cluster_id == 9), 4])
fraq_table["ES55O", "10"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55O" & pyclone_out$cluster_id == 10), 4])
fraq_table["ES55O", "11"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55O" & pyclone_out$cluster_id == 11), 4])
fraq_table["ES55O", "12"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55O" & pyclone_out$cluster_id == 12), 4])
fraq_table["ES55O", "13"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55O" & pyclone_out$cluster_id == 13), 4])
fraq_table["ES55O", "14"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55O" & pyclone_out$cluster_id == 14), 4])
fraq_table["ES55O", "15"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55O" & pyclone_out$cluster_id == 15), 4])
fraq_table["ES55T", "0"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55T" & pyclone_out$cluster_id == 0), 4])
fraq_table["ES55T", "1"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55T" & pyclone_out$cluster_id == 1), 4])
fraq_table["ES55T", "2"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55T" & pyclone_out$cluster_id == 2), 4])
fraq_table["ES55T", "3"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55T" & pyclone_out$cluster_id == 3), 4])
fraq_table["ES55T", "4"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55T" & pyclone_out$cluster_id == 4), 4])
fraq_table["ES55T", "5"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55T" & pyclone_out$cluster_id == 5), 4])
fraq_table["ES55T", "6"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55T" & pyclone_out$cluster_id == 6), 4])
fraq_table["ES55T", "7"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55T" & pyclone_out$cluster_id == 7), 4])
fraq_table["ES55T", "8"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55T" & pyclone_out$cluster_id == 8), 4])
fraq_table["ES55T", "9"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55T" & pyclone_out$cluster_id == 9), 4])
fraq_table["ES55T", "10"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55T" & pyclone_out$cluster_id == 10), 4])
fraq_table["ES55T", "11"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55T" & pyclone_out$cluster_id == 11), 4])
fraq_table["ES55T", "12"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55T" & pyclone_out$cluster_id == 12), 4])
fraq_table["ES55T", "13"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55T" & pyclone_out$cluster_id == 13), 4])
fraq_table["ES55T", "14"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55T" & pyclone_out$cluster_id == 14), 4])
fraq_table["ES55T", "15"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES55T" & pyclone_out$cluster_id == 15), 4])
fraq_table["ES80O", "0"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80O" & pyclone_out$cluster_id == 0), 4])
fraq_table["ES80O", "1"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80O" & pyclone_out$cluster_id == 1), 4])
fraq_table["ES80O", "2"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80O" & pyclone_out$cluster_id == 2), 4])
fraq_table["ES80O", "3"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80O" & pyclone_out$cluster_id == 3), 4])
fraq_table["ES80O", "4"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80O" & pyclone_out$cluster_id == 4), 4])
fraq_table["ES80O", "5"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80O" & pyclone_out$cluster_id == 5), 4])
fraq_table["ES80O", "6"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80O" & pyclone_out$cluster_id == 6), 4])
fraq_table["ES80O", "7"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80O" & pyclone_out$cluster_id == 7), 4])
fraq_table["ES80O", "8"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80O" & pyclone_out$cluster_id == 8), 4])
fraq_table["ES80O", "9"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80O" & pyclone_out$cluster_id == 9), 4])
fraq_table["ES80O", "10"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80O" & pyclone_out$cluster_id == 10), 4])
fraq_table["ES80O", "11"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80O" & pyclone_out$cluster_id == 11), 4])
fraq_table["ES80O", "12"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80O" & pyclone_out$cluster_id == 12), 4])
fraq_table["ES80O", "13"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80O" & pyclone_out$cluster_id == 13), 4])
fraq_table["ES80O", "14"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80O" & pyclone_out$cluster_id == 14), 4])
fraq_table["ES80O", "15"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80O" & pyclone_out$cluster_id == 15), 4])
fraq_table["ES80T", "0"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80T" & pyclone_out$cluster_id == 0), 4])
fraq_table["ES80T", "1"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80T" & pyclone_out$cluster_id == 1), 4])
fraq_table["ES80T", "2"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80T" & pyclone_out$cluster_id == 2), 4])
fraq_table["ES80T", "3"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80T" & pyclone_out$cluster_id == 3), 4])
fraq_table["ES80T", "4"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80T" & pyclone_out$cluster_id == 4), 4])
fraq_table["ES80T", "5"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80T" & pyclone_out$cluster_id == 5), 4])
fraq_table["ES80T", "6"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80T" & pyclone_out$cluster_id == 6), 4])
fraq_table["ES80T", "7"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80T" & pyclone_out$cluster_id == 7), 4])
fraq_table["ES80T", "8"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80T" & pyclone_out$cluster_id == 8), 4])
fraq_table["ES80T", "9"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80T" & pyclone_out$cluster_id == 9), 4])
fraq_table["ES80T", "10"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80T" & pyclone_out$cluster_id == 10), 4])
fraq_table["ES80T", "11"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80T" & pyclone_out$cluster_id == 11), 4])
fraq_table["ES80T", "12"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80T" & pyclone_out$cluster_id == 12), 4])
fraq_table["ES80T", "13"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80T" & pyclone_out$cluster_id == 13), 4])
fraq_table["ES80T", "14"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80T" & pyclone_out$cluster_id == 14), 4])
fraq_table["ES80T", "15"] <- unique(pyclone_out[which(pyclone_out$sample_id == "ES80T" & pyclone_out$cluster_id == 15), 4])

# put the frequency table in the correct format for riverplots
table_river <- as.data.frame(fraq_table)
table_river$cluster <- paste0("cluster_", table_river$cluster_id)
table_river$variable <- table_river$sample_id
table_river$value <- table_river$Freq

############### plot the rivers ####################
ggplot(table_river[which(table_river$variable %in% c("ES16T", "ES55T", "ES80T") & 
                           table_river$cluster %in% c("cluster_1", "cluster_3",
                                        "cluster_4", "cluster_6",
                                        "cluster_7", "cluster_8", "cluster_9",
                                        "cluster_10", "cluster_11")),], 
       aes(x = variable, y = value, stratum = cluster, alluvium = cluster, fill = cluster, label = cluster)) + theme_classic() +
  geom_flow(curve_type = "sigmoid") + scale_x_discrete(expand = c(.1, .1)) + geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3)+ labs(x="Tumor vs 1st relapse vs 2nd relapse", y="cellular prevalence") +
  theme(legend.position = "none") + scale_fill_manual(values = rainbow(9)) + 
  geom_point(aes(y = 1.12, x = 2.55, shape = "1")) + annotate(
    "text", label = "STAG2", x = 2.64, y = 1.14, size = 4, colour = "black") +
  ggtitle("Clonal expansion over different timepoints") + scale_shape_manual(values = c(8))


ggplot(table_river[which(table_river$variable %in% c("ES16T", "ES16O") & 
                           table_river$cluster %in% c("cluster_0", "cluster_1",
                                        "cluster_8", "cluster_9", "cluster_10", "cluster_11",
                                        "cluster_12", "cluster_13")),], 
       aes(x = factor(variable, levels = c("ES16T", "ES16O")), y = value, stratum = cluster, alluvium = cluster, fill = cluster, label = cluster)) + theme_classic() +
  geom_flow(curve_type = "sigmoid") + scale_x_discrete(expand = c(.1, .1)) + geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) + labs(x="Tumor vs Organoid", y="cellular prevalence")+
  theme(legend.position = "none") + scale_fill_manual(values = rainbow(8)) + 
  ggtitle("Clonal expansion in organoid culture for ES-016")

ggplot(table_river[which(table_river$variable %in% c("ES55T", "ES55O") & 
                           table_river$cluster %in% c("cluster_0", "cluster_1", "cluster_3",
                                        "cluster_4", "cluster_5",
                                         "cluster_9", "cluster_10", "cluster_11",
                                        "cluster_12", "cluster_13")),], 
       aes(x = factor(variable, levels = c("ES55T", "ES55O")), y = value, stratum = cluster, alluvium = cluster, fill = cluster, label = cluster)) + theme_classic() +
  geom_flow(curve_type = "sigmoid") + scale_x_discrete(expand = c(.1, .1)) + geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) + labs(x="Tumor vs Organoid", y="cellular prevalence")+
  theme(legend.position = "none") + scale_fill_manual(values = rainbow(10)) + 
  ggtitle("Clonal expansion in organoid culture for ES-055")

ggplot(table_river[which(table_river$variable %in% c("ES80T", "ES80O") & 
                           table_river$cluster %in% c("cluster_0", "cluster_1","cluster_2",
                                        "cluster_6", "cluster_7",
                                         "cluster_9", "cluster_10", "cluster_11",
                                        "cluster_12")),], 
       aes(x = factor(variable, levels = c("ES80T", "ES80O")), y = value, stratum = cluster, alluvium = cluster, fill = cluster, label = cluster)) + theme_classic() +
  geom_flow(curve_type = "sigmoid") + scale_x_discrete(expand = c(.1, .1)) + geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3)+ labs(x="Tumor vs Organoid", y="cellular prevalence") +
  theme(legend.position = "none") + scale_fill_manual(values = rainbow(9)) + 
  ggtitle("Clonal expansion in organoid culture for ES-080")


########## other figure ##########

ggplot(table_river[which(table_river$cluster_id %in% c("0", "1", "2", "3", "4", "5", "6", "7", 
                                                       "8", "9", "10", "11", "12", "13", "15")),]) + 
  geom_col(aes(x = sample_id, y = Freq, fill = cluster_id), position = "fill") + scale_fill_manual(values = rainbow(15)) + theme_light()

picture + theme_minimal() + scale_y_continuous(minor_breaks = c(0:15)) + scale_color_manual(values = rainbow(16))


