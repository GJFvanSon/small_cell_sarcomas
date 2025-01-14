library(maftools)
library(ggplot)
library(gggenes)

maf <- read.maf("fusions_maf.maf", vc_nonSyn = "Gene Fusion", removeDuplicatedVariants = FALSE)
lollipopPlot(maf, gene = "EWSR1", showDomainLabel = F, labelPos = "all", repel = T )
lollipopPlot(maf, gene = "FLI1", showDomainLabel = F, labelPos = "all", repel = T )
lollipopPlot(maf, gene = "CIC", showDomainLabel = F, labelPos = "all", repel = T )
lollipopPlot(maf, gene = "DUX4", showDomainLabel = F, labelPos = "all", repel = T )
lollipopPlot(maf, gene = "CCNB3", showDomainLabel = F, labelPos = "all", repel = F )
lollipopPlot(maf, gene = "BCOR", showDomainLabel = F, labelPos = "all", repel = T )
lollipopPlot(maf, gene = "KMT2D", showDomainLabel = F, labelPos = "all", repel = T )
lollipopPlot(maf, gene = "ERG", showDomainLabel = F, labelPos = "all", repel = T )
lollipopPlot(maf, gene = "FEV", showDomainLabel = F, labelPos = "all", repel = T)



EWSR1_gtf <-read.table("Scripts/EWSR1.gtf", sep = "\t")
colnames(EWSR1_gtf) <- c("Chromosome", "Source", "Type", "Start", "End", "point", "Strand", "something", "Info")
info <- lapply(EWSR1_gtf$Info, function(x){strsplit(x, ";")[[1]]})
EWSR1_gtf$transcript_id <- lapply(info, function(x){if(length(grep("transcript_id", x)) == 0){NA}else{strsplit(x[grep("transcript_id", x)], "id ")[[1]][2]}}) 
EWSR1_gtf$gene_version <- lapply(info, function(x){if(length(grep("gene_version", x)) == 0){NA}else{strsplit(x[grep("gene_version", x)], "on ")[[1]][2]}}) 
EWSR1_gtf$gene_id <- lapply(info, function(x){if(length(grep("gene_id", x)) == 0){NA}else{strsplit(x[grep("gene_id", x)], "id ")[[1]][2]}}) 
EWSR1_gtf$exon_number <- lapply(info, function(x){if(length(grep("exon_number", x)) == 0){NA}else{strsplit(x[grep("exon_number", x)], "er ")[[1]][2]}}) 
EWSR1_gtf$gene_name <- lapply(info, function(x){if(length(grep("gene_name", x)) == 0){NA}else{strsplit(x[grep("gene_name", x)], "me ")[[1]][2]}}) 
unlist(unique(EWSR1_gtf$transcript_id))
EWSR1_gtf_659 <- EWSR1_gtf[c(1,which(EWSR1_gtf$transcript_id == "ENST00000629659")),]

FLI1_gtf <-read.table("Scripts/FLI1.gtf", sep = "\t")
colnames(FLI1_gtf) <- c("Chromosome", "Source", "Type", "Start", "End", "point", "Strand", "something", "Info")
info <- lapply(FLI1_gtf$Info, function(x){strsplit(x, ";")[[1]]})
FLI1_gtf$transcript_id <- lapply(info, function(x){if(length(grep("transcript_id", x)) == 0){NA}else{strsplit(x[grep("transcript_id", x)], "id ")[[1]][2]}}) 
FLI1_gtf$gene_version <- lapply(info, function(x){if(length(grep("gene_version", x)) == 0){NA}else{strsplit(x[grep("gene_version", x)], "on ")[[1]][2]}}) 
FLI1_gtf$gene_id <- lapply(info, function(x){if(length(grep("gene_id", x)) == 0){NA}else{strsplit(x[grep("gene_id", x)], "id ")[[1]][2]}}) 
FLI1_gtf$exon_number <- lapply(info, function(x){if(length(grep("exon_number", x)) == 0){NA}else{strsplit(x[grep("exon_number", x)], "er ")[[1]][2]}}) 
FLI1_gtf$gene_name <- lapply(info, function(x){if(length(grep("gene_name", x)) == 0){NA}else{strsplit(x[grep("gene_name", x)], "me ")[[1]][2]}}) 
unlist(unique(FLI1_gtf$transcript_id))
FLI1_gtf_954 <- FLI1_gtf[c(1,which(FLI1_gtf$transcript_id == "ENST00000344954")),]

ERG_gtf <-read.table("Scripts/ERG.gtf", sep = "\t")
colnames(ERG_gtf) <- c("Chromosome", "Source", "Type", "Start", "End", "point", "Strand", "something", "Info")
info <- lapply(ERG_gtf$Info, function(x){strsplit(x, ";")[[1]]})
ERG_gtf$transcript_id <- lapply(info, function(x){if(length(grep("transcript_id", x)) == 0){NA}else{strsplit(x[grep("transcript_id", x)], "id ")[[1]][2]}}) 
ERG_gtf$gene_version <- lapply(info, function(x){if(length(grep("gene_version", x)) == 0){NA}else{strsplit(x[grep("gene_version", x)], "on ")[[1]][2]}}) 
ERG_gtf$gene_id <- lapply(info, function(x){if(length(grep("gene_id", x)) == 0){NA}else{strsplit(x[grep("gene_id", x)], "id ")[[1]][2]}}) 
ERG_gtf$exon_number <- lapply(info, function(x){if(length(grep("exon_number", x)) == 0){NA}else{strsplit(x[grep("exon_number", x)], "er ")[[1]][2]}}) 
ERG_gtf$gene_name <- lapply(info, function(x){if(length(grep("gene_name", x)) == 0){NA}else{strsplit(x[grep("gene_name", x)], "me ")[[1]][2]}}) 
unlist(unique(ERG_gtf$transcript_id))
ERG_gtf_319 <- ERG_gtf[c(1,which(ERG_gtf$transcript_id == "ENST00000288319")),]

CIC_gtf <-read.table("Scripts/CIC.gtf", sep = "\t")
colnames(CIC_gtf) <- c("Chromosome", "Source", "Type", "Start", "End", "point", "Strand", "something", "Info")
info <- lapply(CIC_gtf$Info, function(x){strsplit(x, ";")[[1]]})
CIC_gtf$transcript_id <- lapply(info, function(x){if(length(grep("transcript_id", x)) == 0){NA}else{strsplit(x[grep("transcript_id", x)], "id ")[[1]][2]}}) 
CIC_gtf$gene_version <- lapply(info, function(x){if(length(grep("gene_version", x)) == 0){NA}else{strsplit(x[grep("gene_version", x)], "on ")[[1]][2]}}) 
CIC_gtf$gene_id <- lapply(info, function(x){if(length(grep("gene_id", x)) == 0){NA}else{strsplit(x[grep("gene_id", x)], "id ")[[1]][2]}}) 
CIC_gtf$exon_number <- lapply(info, function(x){if(length(grep("exon_number", x)) == 0){NA}else{strsplit(x[grep("exon_number", x)], "er ")[[1]][2]}}) 
CIC_gtf$gene_name <- lapply(info, function(x){if(length(grep("gene_name", x)) == 0){NA}else{strsplit(x[grep("gene_name", x)], "me ")[[1]][2]}}) 
unlist(unique(CIC_gtf$transcript_id))
CIC_gtf_681 <- CIC_gtf[c(1,which(CIC_gtf$transcript_id == "ENST00000572681")),]

DUX4_gtf <-read.table("Scripts/DUX4.gtf", sep = "\t")
colnames(DUX4_gtf) <- c("Chromosome", "Source", "Type", "Start", "End", "point", "Strand", "something", "Info")
info <- lapply(DUX4_gtf$Info, function(x){strsplit(x, ";")[[1]]})
DUX4_gtf$transcript_id <- lapply(info, function(x){if(length(grep("transcript_id", x)) == 0){NA}else{strsplit(x[grep("transcript_id", x)], "id ")[[1]][2]}}) 
DUX4_gtf$gene_version <- lapply(info, function(x){if(length(grep("gene_version", x)) == 0){NA}else{strsplit(x[grep("gene_version", x)], "on ")[[1]][2]}}) 
DUX4_gtf$gene_id <- lapply(info, function(x){if(length(grep("gene_id", x)) == 0){NA}else{strsplit(x[grep("gene_id", x)], "id ")[[1]][2]}}) 
DUX4_gtf$exon_number <- lapply(info, function(x){if(length(grep("exon_number", x)) == 0){NA}else{strsplit(x[grep("exon_number", x)], "er ")[[1]][2]}}) 
DUX4_gtf$gene_name <- lapply(info, function(x){if(length(grep("gene_name", x)) == 0){NA}else{strsplit(x[grep("gene_name", x)], "me ")[[1]][2]}}) 
unlist(unique(DUX4_gtf$transcript_id))
DUX4_gtf_263 <- DUX4_gtf[c(1,which(DUX4_gtf$transcript_id == "ENST00000570263")),]

BCOR_gtf <-read.table("Scripts/BCOR.gtf", sep = "\t")
colnames(BCOR_gtf) <- c("Chromosome", "Source", "Type", "Start", "End", "point", "Strand", "something", "Info")
info <- lapply(BCOR_gtf$Info, function(x){strsplit(x, ";")[[1]]})
BCOR_gtf$transcript_id <- lapply(info, function(x){if(length(grep("transcript_id", x)) == 0){NA}else{strsplit(x[grep("transcript_id", x)], "id ")[[1]][2]}}) 
BCOR_gtf$gene_version <- lapply(info, function(x){if(length(grep("gene_version", x)) == 0){NA}else{strsplit(x[grep("gene_version", x)], "on ")[[1]][2]}}) 
BCOR_gtf$gene_id <- lapply(info, function(x){if(length(grep("gene_id", x)) == 0){NA}else{strsplit(x[grep("gene_id", x)], "id ")[[1]][2]}}) 
BCOR_gtf$exon_number <- lapply(info, function(x){if(length(grep("exon_number", x)) == 0){NA}else{strsplit(x[grep("exon_number", x)], "er ")[[1]][2]}}) 
BCOR_gtf$gene_name <- lapply(info, function(x){if(length(grep("gene_name", x)) == 0){NA}else{strsplit(x[grep("gene_name", x)], "me ")[[1]][2]}}) 
unlist(unique(BCOR_gtf$transcript_id))
BCOR_gtf_905 <- BCOR_gtf[c(1,which(BCOR_gtf$transcript_id == "ENST00000413905")),]

KMT2D_gtf <-read.table("Scripts/KMT2D.gtf", sep = "\t")
colnames(KMT2D_gtf) <- c("Chromosome", "Source", "Type", "Start", "End", "point", "Strand", "something", "Info")
info <- lapply(KMT2D_gtf$Info, function(x){strsplit(x, ";")[[1]]})
KMT2D_gtf$transcript_id <- lapply(info, function(x){if(length(grep("transcript_id", x)) == 0){NA}else{strsplit(x[grep("transcript_id", x)], "id ")[[1]][2]}}) 
KMT2D_gtf$gene_version <- lapply(info, function(x){if(length(grep("gene_version", x)) == 0){NA}else{strsplit(x[grep("gene_version", x)], "on ")[[1]][2]}}) 
KMT2D_gtf$gene_id <- lapply(info, function(x){if(length(grep("gene_id", x)) == 0){NA}else{strsplit(x[grep("gene_id", x)], "id ")[[1]][2]}}) 
KMT2D_gtf$exon_number <- lapply(info, function(x){if(length(grep("exon_number", x)) == 0){NA}else{strsplit(x[grep("exon_number", x)], "er ")[[1]][2]}}) 
KMT2D_gtf$gene_name <- lapply(info, function(x){if(length(grep("gene_name", x)) == 0){NA}else{strsplit(x[grep("gene_name", x)], "me ")[[1]][2]}}) 
unlist(unique(KMT2D_gtf$transcript_id))
KMT2D_gtf_067 <- KMT2D_gtf[c(1,which(KMT2D_gtf$transcript_id == "ENST00000301067")),]

CCNB3_gtf <-read.table("Scripts/CCNB3.gtf", sep = "\t")
colnames(CCNB3_gtf) <- c("Chromosome", "Source", "Type", "Start", "End", "point", "Strand", "something", "Info")
info <- lapply(CCNB3_gtf$Info, function(x){strsplit(x, ";")[[1]]})
CCNB3_gtf$transcript_id <- lapply(info, function(x){if(length(grep("transcript_id", x)) == 0){NA}else{strsplit(x[grep("transcript_id", x)], "id ")[[1]][2]}}) 
CCNB3_gtf$gene_version <- lapply(info, function(x){if(length(grep("gene_version", x)) == 0){NA}else{strsplit(x[grep("gene_version", x)], "on ")[[1]][2]}}) 
CCNB3_gtf$gene_id <- lapply(info, function(x){if(length(grep("gene_id", x)) == 0){NA}else{strsplit(x[grep("gene_id", x)], "id ")[[1]][2]}}) 
CCNB3_gtf$exon_number <- lapply(info, function(x){if(length(grep("exon_number", x)) == 0){NA}else{strsplit(x[grep("exon_number", x)], "er ")[[1]][2]}}) 
CCNB3_gtf$gene_name <- lapply(info, function(x){if(length(grep("gene_name", x)) == 0){NA}else{strsplit(x[grep("gene_name", x)], "me ")[[1]][2]}}) 
unlist(unique(CCNB3_gtf$transcript_id))
CCNB3_gtf_042 <- CCNB3_gtf[c(1,which(CCNB3_gtf$transcript_id == "ENST00000376042")),]

FEV_gtf <-read.table("Scripts/FEV.gtf", sep = "\t")
colnames(FEV_gtf) <- c("Chromosome", "Source", "Type", "Start", "End", "point", "Strand", "something", "Info")
info <- lapply(FEV_gtf$Info, function(x){strsplit(x, ";")[[1]]})
FEV_gtf$transcript_id <- lapply(info, function(x){if(length(grep("transcript_id", x)) == 0){NA}else{strsplit(x[grep("transcript_id", x)], "id ")[[1]][2]}}) 
FEV_gtf$gene_version <- lapply(info, function(x){if(length(grep("gene_version", x)) == 0){NA}else{strsplit(x[grep("gene_version", x)], "on ")[[1]][2]}}) 
FEV_gtf$gene_id <- lapply(info, function(x){if(length(grep("gene_id", x)) == 0){NA}else{strsplit(x[grep("gene_id", x)], "id ")[[1]][2]}}) 
FEV_gtf$exon_number <- lapply(info, function(x){if(length(grep("exon_number", x)) == 0){NA}else{strsplit(x[grep("exon_number", x)], "er ")[[1]][2]}}) 
FEV_gtf$gene_name <- lapply(info, function(x){if(length(grep("gene_name", x)) == 0){NA}else{strsplit(x[grep("gene_name", x)], "me ")[[1]][2]}}) 
unlist(unique(FEV_gtf$transcript_id))
FEV_gtf_521 <- FEV_gtf[c(1,which(FEV_gtf$transcript_id == "rna-NM_017521.3.1")),]


for_figure <- data.frame("molecule" = unlist(c(EWSR1_gtf_659[which(EWSR1_gtf_659$Type == "exon"),]$gene_name,
                                               CIC_gtf_681[which(CIC_gtf_681$Type == "exon"),]$gene_name,
                                               DUX4_gtf_263[which(DUX4_gtf_263$Type == "exon"),]$gene_name,
                                               BCOR_gtf_905[which(BCOR_gtf_905$Type == "exon"),]$gene_name,
                                               KMT2D_gtf_067[which(KMT2D_gtf_067$Type == "exon"),]$gene_name,
                                               CCNB3_gtf_042[which(CCNB3_gtf_042$Type == "exon"),]$gene_name,
                                               ERG_gtf_319[which(ERG_gtf_319$Type == "exon"),]$gene_name,
                                               FLI1_gtf_954[which(FLI1_gtf_954$Type == "exon"),]$gene_name,
                                               FEV_gtf_521[which(FEV_gtf_521$Type == "exon"),]$gene_name)),
                         "exon" = unlist(c(EWSR1_gtf_659[which(EWSR1_gtf_659$Type == "exon"),]$exon_number,
                                           CIC_gtf_681[which(CIC_gtf_681$Type == "exon"),]$exon_number,
                                           DUX4_gtf_263[which(DUX4_gtf_263$Type == "exon"),]$exon_number,
                                           BCOR_gtf_905[which(BCOR_gtf_905$Type == "exon"),]$exon_number,
                                           KMT2D_gtf_067[which(KMT2D_gtf_067$Type == "exon"),]$exon_number,
                                           CCNB3_gtf_042[which(CCNB3_gtf_042$Type == "exon"),]$exon_number,
                                           ERG_gtf_319[which(ERG_gtf_319$Type == "exon"),]$exon_number,
                                           FLI1_gtf_954[which(FLI1_gtf_954$Type == "exon"),]$exon_number,
                                           FEV_gtf_521[which(FEV_gtf_521$Type == "exon"),]$exon_number)),
                         "start" = unlist(c(EWSR1_gtf_659[which(EWSR1_gtf_659$Type == "exon"),]$Start,
                                            CIC_gtf_681[which(CIC_gtf_681$Type == "exon"),]$Start,
                                            DUX4_gtf_263[which(DUX4_gtf_263$Type == "exon"),]$Start,
                                            BCOR_gtf_905[which(BCOR_gtf_905$Type == "exon"),]$Start,
                                            KMT2D_gtf_067[which(KMT2D_gtf_067$Type == "exon"),]$Start,
                                            CCNB3_gtf_042[which(CCNB3_gtf_042$Type == "exon"),]$Start,
                                            ERG_gtf_319[which(ERG_gtf_319$Type == "exon"),]$Start,
                                            FLI1_gtf_954[which(FLI1_gtf_954$Type == "exon"),]$Start,
                                            FEV_gtf_521[which(FEV_gtf_521$Type == "exon"),]$Start)),
                         "end" =unlist(c(EWSR1_gtf_659[which(EWSR1_gtf_659$Type == "exon"),]$End,
                                         CIC_gtf_681[which(CIC_gtf_681$Type == "exon"),]$End,
                                         DUX4_gtf_263[which(DUX4_gtf_263$Type == "exon"),]$End,
                                         BCOR_gtf_905[which(BCOR_gtf_905$Type == "exon"),]$End,
                                         KMT2D_gtf_067[which(KMT2D_gtf_067$Type == "exon"),]$End,
                                         CCNB3_gtf_042[which(CCNB3_gtf_042$Type == "exon"),]$End,
                                         ERG_gtf_319[which(ERG_gtf_319$Type == "exon"),]$End,
                                         FLI1_gtf_954[which(FLI1_gtf_954$Type == "exon"),]$End,
                                         FEV_gtf_521[which(FEV_gtf_521$Type == "exon"),]$End)),
                         "strand" = unlist(c(EWSR1_gtf_659[which(EWSR1_gtf_659$Type == "exon"),]$Strand,
                                             CIC_gtf_681[which(CIC_gtf_681$Type == "exon"),]$Strand,
                                             DUX4_gtf_263[which(DUX4_gtf_263$Type == "exon"),]$Strand,
                                             BCOR_gtf_905[which(BCOR_gtf_905$Type == "exon"),]$Strand,
                                             KMT2D_gtf_067[which(KMT2D_gtf_067$Type == "exon"),]$Strand,
                                             CCNB3_gtf_042[which(CCNB3_gtf_042$Type == "exon"),]$Strand,
                                             ERG_gtf_319[which(ERG_gtf_319$Type == "exon"),]$Strand,
                                             FLI1_gtf_954[which(FLI1_gtf_954$Type == "exon"),]$Strand,
                                             FEV_gtf_521[which(FEV_gtf_521$Type == "exon"),]$Strand)),
                         "orientation" = 0)


feature_data <- data.frame("position" = maf@data$Start_Position, "molecule" = maf@data$Hugo_Symbol, 
                           "name" = maf@data$Tumor_Sample_Barcode, "forward" = T)
ggplot(for_figure, aes(xmin = start, xmax = end, y = molecule, fill = "Exon")) +
  geom_feature(
    data = feature_data,
    aes(x = position, y = molecule, forward = forward), feature_width = unit(-2, "mm")
  ) + scale_fill_manual(values = c("steelblue3")) +
  geom_gene_arrow(arrow_body_height = unit(5, "mm"), arrowhead_height = unit(5, "mm"), arrowhead_width = unit(0, "mm")) +
  geom_label_repel(data = feature_data, aes(x = position, label = name, y = molecule), size = 2.6, inherit.aes = F,
                   arrow = arrow(length = unit(0.02, "npc")), nudge_y = -1) + 
  facet_wrap(~ factor(molecule, levels = c("ERG" , "EWSR1", "FLI1", "CIC", "DUX4", "CCNB3", "BCOR", "KMT2D", "FEV")), scales = "free", ncol = 1) + 
  theme_genes() +NoLegend() + ggtitle("Breakpoints in fusion genes")
