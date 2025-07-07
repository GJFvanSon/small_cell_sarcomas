############################################################################################
#                  RNA bulk analysis with chip-seq on Ewing(-like) sarcomas                  #
#                                                                                          #
#                                   Gijs van Son                                           #          
#                                   June 11 2024                                           #
#                         g.j.f.vanson@prinsesmaximacentrum.nl                             #                            
#                                                                                          #
############################################################################################
########### setting the environment #############
library(data.table)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggvenn)
library(enrichR)
library(Seurat)
library(readxl)
library(viridis)
setwd("/Users/g.j.f.vanson/OneDrive - Prinses Maxima Centrum/Femke/CIC_DUX4/RNA_BULK/")

############ count data files processing ####################
# files <- list.files("New_count_files/", full.names = T)
# 
# 
# 
# file_list <- lapply(files, function(x){fread(x)})
# file_list[[1]]
# count_data <- data.frame("GeneID" = file_list[[1]]$GeneName)
# for(x in c(1:length(file_list))){
#   name <- strsplit(strsplit(files[x], "//")[[1]][2],"_")[[1]][1]
#   message(name)
#   count_data[[name]] <- file_list[[x]][[2]]
# }
# count_data
# count_data_dup <- count_data[which(duplicated(count_data$GeneID)),]
# for(gene in unique(count_data_dup$GeneID)){
#   count_data <- count_data[-which(count_data$GeneID == gene)[-1],]
#   count_data[which(count_data$GeneID == gene),-1] <- 
#     colSums(count_data_dup[which(count_data_dup$GeneID == gene),-1])
# }
# rownames(count_data) <- count_data$GeneID
# count_data <- count_data[,-1]
# write.table(count_data, "readcounts_table.csv")
######### preparing and running DESEQ2 ####################
#input data
count_data <- read.table("readcounts_table.csv") # save this in work dir
#file_list <- NULL
#count_data_dup <- NULL
count_data <- count_data[,-grep("PMOBM000APM", colnames(count_data))]
metadata <- read.table("metadata2.csv", sep = ",", header = T, row.names = 1) # save this in workdir
# checks and prep for deseq
all(rownames(metadata) %in% colnames(count_data))
all(colnames(count_data) %in% rownames(metadata))
meta_data <- metadata[colnames(count_data),]
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = meta_data, design = ~ 0+Fusion+Type)
featureData <- data.frame(gene=rownames(count_data))
mcols(dds) <- DataFrame(mcols(dds), featureData)
dds$Fusion <- factor(dds$Fusion, levels = c("EWS-FLI1", "RMS", "BCOR","EWS-ERG","EWS-FEV", "CIC-DUX4"))

# running deseq
dds <- DESeq(dds)

res <- results(dds, contrast = c("Fusion", "EWS-FLI1", "CIC-DUX4"))
ntd <- normTransform(dds)
normalised_counts <- assay(ntd)
vsd <- vst(dds, blind=FALSE)
res_significant <- res[which(res$pvalue < 0.05),]
interesting_rows1 <- unlist(list(order(res_significant$log2FoldChange, decreasing =TRUE)[1:100]))
interesting_rows2 <-  unlist(list(order(res_significant$log2FoldChange, decreasing =F)[1:100]))

######## plotting DESEQ2 results ############
heatmap <- pheatmap(assay(ntd)[c(rownames(res_significant)[interesting_rows1][c(1:20)], rownames(res_significant)[interesting_rows2][c(1:20)]),], 
                    cluster_rows = F, 
                    show_rownames = T, 
                    cluster_cols = T, 
                    color = colorRampPalette(rev(brewer.pal(8, name = "PuOr")))(255),
                    border_color = "white", 
                    scale = "row",
                    main = "RNA expression data of differentially expressed genes of EWS-FLI vs CIC-DUX4",
                    annotation_col = meta_data[,c("Fusion", "Patient_ID")])
sampleDist <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDist)
colors <- colorRampPalette(rev(brewer.pal(11, name = "PuOr")))(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDist, 
         clustering_distance_cols=sampleDist,
         annotation_col = meta_data[,c("Fusion", "Patient_ID")],
         col=viridis(100),
         main = "sample distance matrix",
         treeheight_col = 0)
# select which samples to show
ordercols <- colnames(assay(ntd))[c(grep("EWS-FLI", meta_data$Fusion),grep("EWS-ERG", meta_data$Fusion),
                                    grep("EWS-FEV", meta_data$Fusion),grep("BCOR", meta_data$Fusion),
                                    grep("CIC-DUX4", meta_data$Fusion),grep("RMS", meta_data$Fusion))]
# remove the RMS samples from plots and the ES10 and ES84 samples that are not so OK
ordercols <- ordercols[-c(15, 54)]
ordercols <- ordercols[-c(54:59)]
pheatmap(assay(vsd)[c("CD99", "CCND1", "IGF1", "PAX7", "NEUROD6", "NKX2-2", "JAK1", "CAV1", "HOXD13", "KDSR", "BCL11B", "GLG1", "ATP1A1", "XAGE2",
                      "ERG", "VIM", "DAD1", "XRCC5", "PRKDC", "BCOR", "BCORL1", "CCNB3", "TLE1", "SATB2", "IGF2",
                      "BCL2", "NTRK3", "ROBO1", "CCND2", "MAML2", "VCAN","FEV" ,"BCL6", "GREM1", "IGF1R", "COL6A1", "JAK2", 
                      "WT1", "CIC", "DUSP4", "ETS1", "ETV1", "ETV4", "ETV5", "DUX4", "MCL1", "MUC5AC", "BCL2", "DUSP6", "HMGA2", "SPRED1",
                      "SPRY4", "VGF", "COL6A5", "COL6A2", "TGFBR2"),
                    ordercols],
         cluster_rows = F,
         show_rownames = T,
         cluster_cols = F,
         color = viridis(5),
         scale = "row",
         border_color = "gray80", 
         annotation_col = meta_data[, c("Fusion", "Type")],
         labels_col = meta_data[ordercols,5])
# pca plotje op twee manieren
plotPCA(vsd, intgroup = c("Fusion"), ntop = 5000) + theme_bw() + 
  scale_color_manual(values = rainbow(6)) + theme_classic() + 
  ggtitle("PCA plot of samples ")
plotPCA(vsd, intgroup = c("Type"), ntop = 5000) + theme_bw() + 
  scale_color_manual(values = rainbow(3)) + theme_classic() + 
  ggtitle("PCA plot of samples ")
PCA_data <- plotPCA(vsd, intgroup = c("Fusion", "Type", "ES_number"), ntop = 5000)
PCA_data <- PCA_data$data[-grep("RMS", PCA_data$data$Fusion),]
PCA_data <- PCA_data[-grep("delete", PCA_data$ES_number),]
ggplot(PCA_data) + geom_point(aes(x = PC1, y=PC2, color=Fusion, shape = Type, size = "test")) + theme_light() +
  scale_size_manual(values = c(3)) + scale_color_manual(values = c("firebrick", "midnightblue", "seagreen",
                                                                   "coral2", "yellow3")) +
  ggtitle("PCA plot of samples") + labs(x = "PC1: 24% variance", y="PC2: 13% variance")
########## chip seq processing ###########
#chipseqdata erbij halen
# https://www.pnas.org/doi/10.1073/pnas.2009137117#supplementary-materials 
# within 3Kb of TSS
# CIC_dux_peaks
chip_peaks_cic <- fread("../chip_seq/chip_seq_peaks_cic_dux4.tsv")
unique(chip_peaks_cic$SYMBOL)[which(unique(chip_peaks_cic$SYMBOL) %in% rownames(res_significant))]

#EWS-FLI peaks
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4492343/#SD3
chip_peaks_EWS <- fread("../chip_seq/chip_seq_peaks_ews_fli1.tsv")

#EWS-ERG peaks
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2905244/#S1 
chip_peaks_ERG <- read_xls("../chip_seq/EWS_ERG_peaks.xls")

#BCOR peaks
chip_peak_BCOR <- fread("../chip_seq/New_analysis/BCOR_ENCODE.bed")
colnames(chip_peak_BCOR) <- c("Chrom", "Start", "End", "peak_id", "Signal")

new_table <- data.frame("GeneID" = file_list[[1]]$GeneName, "Chromosome" = file_list[[1]]$Chr, "start" = file_list[[1]]$Start)
new_table$start[grep(";", new_table$start)] <- unlist(lapply(new_table$start[grep(";", new_table$start)], function(x){strsplit(x,";")[[1]][1]}))
new_table$Chromosome[grep(";", new_table$Chromosome)] <- unlist(lapply(new_table$Chromosome[grep(";", new_table$Chromosome)], function(x){strsplit(x,";")[[1]][1]}))
new_table$peak_ID <- NA
new_table$distance <- NA
test_erg <- data.frame("peak_id" = c(1:length(chip_peaks_ERG$Chromosome)),"Chrom" = chip_peaks_ERG$Chromosome, "Start" =chip_peaks_ERG$Start, "End" =chip_peaks_ERG$End, 
                       "signal_FE_overIgg" = chip_peaks_ERG$Fold_change, "P_value" = chip_peaks_ERG$`p-value`, 
                       "Gene" = NA, "GeneTSS" = NA, "GeneChr" = NA, "distance" = NA)
id <- test_erg$peak_id[1]
for(id in test_erg$peak_id){
  id_chrom = paste0("chr",test_erg$Chrom[id])
  id_start = test_erg$Start[id]
  id_end = test_erg$End[id]
  search_table <- new_table[new_table$Chromosome == id_chrom,]
  search_table2 <- search_table[intersect(which(as.numeric(search_table$start) > (id_start - 3000)), which(as.numeric(search_table$start) < (id_end + 3000))),]
  for(gene in search_table2$GeneID){
    new_table$distance[which(new_table$GeneID == gene)] <- min(abs(id_start - as.numeric(search_table2$start[which(search_table2$GeneID == gene)])),
                                                               abs(id_end - as.numeric(search_table2$start[which(search_table2$GeneID == gene)])))
    new_table$peak_ID[which(new_table$GeneID == gene)] <- id
  }
}
EWSERG_chip <- new_table[-which(is.na(new_table$peak_ID)),]
EWSERG_genes <- unique(EWSERG_chip$GeneID)
new_table <- data.frame("GeneID" = file_list[[1]]$GeneName, "Chromosome" = file_list[[1]]$Chr, "start" = file_list[[1]]$Start)
new_table$start[grep(";", new_table$start)] <- unlist(lapply(new_table$start[grep(";", new_table$start)], function(x){strsplit(x,";")[[1]][1]}))
new_table$Chromosome[grep(";", new_table$Chromosome)] <- unlist(lapply(new_table$Chromosome[grep(";", new_table$Chromosome)], function(x){strsplit(x,";")[[1]][1]}))
new_table$peak_ID <- NA
new_table$distance <- NA
for(row in c(1:nrow(chip_peak_BCOR))){
  id_chrom = paste0("chr",chip_peak_BCOR$Chrom[row])
  id_start = chip_peak_BCOR$Start[row]
  id_end = chip_peak_BCOR$End[row]
  search_table <- new_table[new_table$Chromosome == id_chrom,]
  search_table2 <- search_table[intersect(which(as.numeric(search_table$start) > (id_start - 3000)), which(as.numeric(search_table$start) < (id_end + 3000))),]
  for(gene in search_table2$GeneID){
    new_table$distance[which(new_table$GeneID == gene)] <- min(abs(id_start - as.numeric(search_table2$start[which(search_table2$GeneID == gene)])),
                                                               abs(id_end - as.numeric(search_table2$start[which(search_table2$GeneID == gene)])))
    new_table$peak_ID[which(new_table$GeneID == gene)] <- chip_peak_BCOR$peak_id[row]
  }
}
# collect genenames from chipseq datasets
BCOR_chip <- new_table[-which(is.na(new_table$peak_ID)),]
BCOR_genes <- unique(BCOR_chip$GeneID)
CICDUX_genes <- unique(chip_peaks_cic$SYMBOL)
EWSFLI_genes <- unique(chip_peaks_EWS$Gene_id)

# overlap gene names from chipseq whith differentially expressed genes.
res_ewing_fli <- results(dds, name = "FusionEWS.FLI1")
res_significant <- res_ewing_fli[which(res_ewing_fli$padj < 0.01),]
ewing_fli_sign_res <- res_significant
interesting_genes_EWSFLI <- res_significant[intersect(rownames(res_significant), EWSFLI_genes),]
interesting_genes_EWSFLI <- interesting_genes_EWSFLI[which(interesting_genes_EWSFLI$log2FoldChange > 0),]
interesting_genes_EWSFLI <- interesting_genes_EWSFLI[order(interesting_genes_EWSFLI$padj),]

res_ewing_cic <- results(dds, name = "FusionCIC.DUX4")
res_significant <- res_ewing_cic[which(res_ewing_cic$padj < 0.05),]
cic_dux_sign_res <- res_significant
interesting_genes_CICDUX <- res_significant[intersect(rownames(res_significant), CICDUX_genes),]
interesting_genes_CICDUX <- interesting_genes_CICDUX[which(interesting_genes_CICDUX$log2FoldChange > 0),]
interesting_genes_CICDUX <- interesting_genes_CICDUX[order(interesting_genes_CICDUX$padj),]

res_ewing_erg <- results(dds, name = "FusionEWS.ERG")
res_significant <- res_ewing_erg[which(res_ewing_erg$padj < 0.05),]
ewing_erg_sign_res <- res_significant
interesting_genes_EWSERG <- res_significant[intersect(rownames(res_significant), EWSERG_genes),]
interesting_genes_EWSERG <- interesting_genes_EWSERG[which(interesting_genes_EWSERG$log2FoldChange > 0),]
interesting_genes_EWSERG <- interesting_genes_EWSERG[order(interesting_genes_EWSERG$log2FoldChange, decreasing = T),]

res_ewing_bcor <- results(dds, name = "FusionBCOR")
res_significant <- res_ewing_bcor[which(res_ewing_bcor$padj < 0.05),]
bcor_sign_res <- res_significant
interesting_genes_BCOR <- res_significant[intersect(rownames(res_significant), BCOR_genes),]
interesting_genes_BCOR <- interesting_genes_BCOR[which(interesting_genes_BCOR$log2FoldChange > 0),]
interesting_genes_BCOR <- interesting_genes_BCOR[order(interesting_genes_BCOR$padj),]
######### chip seq plots ##############

# this is the plot you have

pheatmap(assay(ntd)[c(rownames(interesting_genes_EWSFLI)[c(1:20)],
                      rownames(interesting_genes_EWSERG)[c(1:20)],
                      rownames(interesting_genes_BCOR)[c(1:20)],
                      rownames(interesting_genes_CICDUX)[c(1:20)], "ETS1", "ETV5", "DUSP4", "IL31RA"),ordercols],
         cluster_rows = F, 
         show_rownames = T, 
         cluster_cols = F, 
         color = viridis(5),
         border_color = "gray90", 
         scale = "row",
         labels_col = meta_data[ordercols,5],
         main = "RNA expression data of differentially expressed genes at chip-seq peaks",
         annotation_col = meta_data[,c("Fusion", "Type")])

############## GOterm analysis ################
# check out go analysis for the same geneset
enrichR::listEnrichrDbs()$libraryName
EWS_FLI_enrich <- as.data.frame(enrichr(rownames(interesting_genes_EWSFLI)[c(1:100)], "GO_Molecular_Function_2023"))[c(1:8),]
EWS_FLI_enrich$fusion <- "EWSR1-FLI1"
EWS_ERG_enrich <- as.data.frame(enrichr(rownames(interesting_genes_EWSERG)[c(1:100)], "GO_Molecular_Function_2023"))[c(99,75),]
EWS_ERG_enrich$fusion <- "EWSR1-ERG"
CIC_DUX_enrich <- as.data.frame(enrichr(rownames(interesting_genes_CICDUX)[c(1:100)], "GO_Molecular_Function_2023"))[c(1:8),]
CIC_DUX_enrich$fusion <- "CIC-DUX4"
BCOR_enrich <- as.data.frame(enrichr(rownames(interesting_genes_BCOR)[c(1:100)], "GO_Molecular_Function_2023"))[c(1:4),]
BCOR_enrich$fusion <- "BCOR"
DE_fusions <- rbind(EWS_FLI_enrich, EWS_ERG_enrich)
DE_fusions <- rbind(DE_fusions, CIC_DUX_enrich)
DE_fusions <- rbind(DE_fusions, BCOR_enrich)
colnames(DE_fusions)
ggplot(DE_fusions) + geom_col(aes(x = factor(GO_Molecular_Function_2023.Term, levels = DE_fusions$GO_Molecular_Function_2023.Term ), y = -log(GO_Molecular_Function_2023.Adjusted.P.value),fill = fusion )) + 
  coord_flip() + theme_classic()
ggplot(DE_fusions) + geom_point(aes(x = fusion, y = GO_Molecular_Function_2023.Term, color = -log(GO_Molecular_Function_2023.Adjusted.P.value), size = -log(GO_Molecular_Function_2023.Adjusted.P.value))) + 
  theme_classic() + scale_color_viridis_c() +  grids(axis = c("xy"), color = "grey92", size = NULL, linetype = NULL) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# check out go analysis for EWS-FLI vs CIC-DUX and EWS-FLI vs EWS-ERG and EWS-FLI vs BCOR
ewsfli_vs_cicdux <- results(dds, contrast = c("Fusion", "EWS.FLI1", "CIC.DUX4"))
ewsfli_vs_cicdux_sign <- ewsfli_vs_cicdux[which(ewsfli_vs_cicdux$padj < 0.05),]
ewsfli_vs_cicdux_sign_pos <- ewsfli_vs_cicdux_sign[which(ewsfli_vs_cicdux_sign$log2FoldChange > 0),]
ewsfli_vs_cicdux_sign_neg <- ewsfli_vs_cicdux_sign[which(ewsfli_vs_cicdux_sign$log2FoldChange < 0),]

ewsfli_vs_ewserg <- results(dds, contrast = c("Fusion", "EWS.FLI1", "EWS.ERG"))
ewsfli_vs_ewserg_sign <- ewsfli_vs_ewserg[which(ewsfli_vs_ewserg$padj < 0.05),]
ewsfli_vs_ewserg_sign_pos <- ewsfli_vs_ewserg_sign[which(ewsfli_vs_ewserg_sign$log2FoldChange > 0),]
ewsfli_vs_ewserg_sign_neg <- ewsfli_vs_ewserg_sign[which(ewsfli_vs_ewserg_sign$log2FoldChange < 0),]

ewsfli_vs_bcor <- results(dds, contrast = c("Fusion", "EWS.FLI1", "BCOR"))
ewsfli_vs_bcor_sign <- ewsfli_vs_cicdux[which(ewsfli_vs_bcor$padj < 0.05),]
ewsfli_vs_bcor_sign_pos <- ewsfli_vs_cicdux_sign[which(ewsfli_vs_bcor_sign$log2FoldChange > 0),]
ewsfli_vs_bcor_sign_neg <- ewsfli_vs_cicdux_sign[which(ewsfli_vs_bcor_sign$log2FoldChange < 0),]

CIC_enrich_up <- as.data.frame(enrichr(rownames(ewsfli_vs_cicdux_sign_pos[order(ewsfli_vs_cicdux_sign_pos$padj)[c(1:100)],]), 
                                       "GO_Molecular_Function_2023"))[c(1:10),]
CIC_enrich_up$fusion <- "CIC-DUX_upregulated"
CIC_enrich_down <- as.data.frame(enrichr(rownames(ewsfli_vs_cicdux_sign_neg[order(ewsfli_vs_cicdux_sign_neg$padj)[c(1:100)],]), 
                                       "GO_Molecular_Function_2023"))[c(1:10),]
CIC_enrich_down$fusion <- "CIC-DUX_downregulated"
ERG_enrich_up <- as.data.frame(enrichr(rownames(ewsfli_vs_ewserg_sign_pos[order(ewsfli_vs_ewserg_sign_pos$padj)[c(1:100)],]), 
                                       "GO_Molecular_Function_2023"))[c(1:10),]
ERG_enrich_up$fusion <- "EWSR1-ERG_upregulated"
ERG_enrich_down <- as.data.frame(enrichr(rownames(ewsfli_vs_ewserg_sign_neg[order(ewsfli_vs_ewserg_sign_neg$padj)[c(1:100)],]), 
                                       "GO_Molecular_Function_2023"))[c(1:10),]
ERG_enrich_down$fusion <- "EWSR1-ERG_downregulated"
BCOR_enrich_up <- as.data.frame(enrichr(rownames(ewsfli_vs_bcor_sign_pos[order(ewsfli_vs_bcor_sign_pos$padj)[c(1:100)],]), 
                                       "GO_Molecular_Function_2023"))[c(1:10),]
BCOR_enrich_up$fusion <- "BCOR_upregulated"
BCOR_enrich_down <- as.data.frame(enrichr(rownames(ewsfli_vs_bcor_sign_pos[order(ewsfli_vs_bcor_sign_pos$padj)[c(1:100)],]), 
                                       "GO_Molecular_Function_2023"))[c(1:10),]
BCOR_enrich_down$fusion <- "BCOR_downregulated"

fusions_up <- rbind(CIC_enrich_up, ERG_enrich_up)
fusions_up <- rbind(fusions_up, BCOR_enrich_up)
  
fusions_down <- rbind(CIC_enrich_down, ERG_enrich_down)
fusions_down <- rbind(fusions_down, BCOR_enrich_down)
ggplot(fusions_up) + geom_point(aes(x = fusion, y = GO_Molecular_Function_2023.Term, color = -log(GO_Molecular_Function_2023.P.value), size = -log(GO_Molecular_Function_2023.P.value))) + 
  theme_classic() + scale_color_viridis_c() +  grids(axis = c("xy"), color = "grey92", size = NULL, linetype = NULL) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(fusions_down) + geom_point(aes(x = fusion, y = GO_Molecular_Function_2023.Term, color = -log(GO_Molecular_Function_2023.P.value), size = -log(GO_Molecular_Function_2023.P.value))) + 
  theme_classic() + scale_color_viridis_c() +  grids(axis = c("xy"), color = "grey92", size = NULL, linetype = NULL) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

########### Venn diagrams ############
# make a venn diagram
genes_ewing_fli <- rownames(ewing_fli_sign_res[which(ewing_fli_sign_res$log2FoldChange > 0),])
genes_ewing_erg <- rownames(ewing_erg_sign_res[which(ewing_erg_sign_res$log2FoldChange > 0),])
genes_cic_dux <- rownames(cic_dux_sign_res[which(cic_dux_sign_res$log2FoldChange > 0),])
genes_bcor <- rownames(bcor_sign_res[which(bcor_sign_res$log2FoldChange > 0),])

genes <- list(genes_ewing_fli, genes_ewing_erg, genes_cic_dux, genes_bcor)
names(genes) <- c("EWSR1-FLI1", "EWSR1-ERG", "CIC-DUX4", "BCOR")

ggvenn(genes, fill_color = rainbow(4)) + ggtitle("significantly upregulated genes vs RMS tumors")

genes_ewing_fli <- rownames(ewing_fli_sign_res[which(ewing_fli_sign_res$log2FoldChange < 0),])
genes_ewing_erg <- rownames(ewing_erg_sign_res[which(ewing_erg_sign_res$log2FoldChange < 0),])
genes_cic_dux <- rownames(cic_dux_sign_res[which(cic_dux_sign_res$log2FoldChange < 0),])
genes_bcor <- rownames(bcor_sign_res[which(bcor_sign_res$log2FoldChange < 0),])

genes <- list(genes_ewing_fli, genes_ewing_erg, genes_cic_dux, genes_bcor)
names(genes) <- c("EWSR1-FLI1", "EWSR1-ERG", "CIC-DUX4", "BCOR")

ggvenn(genes, fill_color = rainbow(4)) + ggtitle("significantly downregulated genes vs RMS tumors")

# have a look at the gene tables in R"
View(as.data.frame(ewing_fli_sign_res))
View(as.data.frame(ewing_erg_sign_res))
View(as.data.frame(cic_dux_sign_res))
View(as.data.frame(bcor_sign_res))









