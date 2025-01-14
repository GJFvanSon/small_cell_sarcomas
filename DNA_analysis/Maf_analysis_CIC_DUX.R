############################################################################################
#                  DNA mutation analysis on Ewing(-like) sarcomas                          #
#                                                                                          #
#                                   Gijs van Son                                           #          
#                                   June 11 2024                                           #
#                         g.j.f.vanson@prinsesmaximacentrum.nl                             #                            
#                                                                                          #
############################################################################################
########### setting the environment #############
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(maftools)
library(readxl)
library(NMF)
library(viridis)
library(gggenes)
library(ggrepel)
library(Seurat)
library(pheatmap)
# Windows
#setwd("D://OneDrive - Prinses Maxima Centrum/Femke/CIC_DUX4/")
# Macbook
setwd("/Users/gijsvanson/OneDrive - Prinses Maxima Centrum/Femke/CIC_DUX4/")

######### MAF-tools analysis #######################
# select variants that are considered synonymous
syn <- c("synonymous_variant", "intron_variant", "", "inframe_deletion_insertion",
         "sequence_feature", "non_coding_transcript_exon_variant", "UTR_variant", "TFBS_variant")
# read the maf-file as a dataframe to get some info from it
df <- data.table::fread("DNA_analysis/MAF_files_all_sarcomas/all_annotated.maf")
#fwrite(df, "DNA_analysis/MAF_files_all_sarcomas/all_annotated.maf", sep = "\t")
# get all the different variant classifications from the file
vc <- c(names(table(df$Variant_Classification)), "gene_fusion", "loss", "gain", "promoter_variant")
# set some variables for reading the maf-file
nonSyn <- setdiff(vc,syn)
colors <- rainbow(length(nonSyn))
names(colors) <- nonSyn
#reading the maffile 
maf <- read.maf("DNA_analysis/MAF_files_all_sarcomas/all_annotated.maf", vc_nonSyn = nonSyn, removeDuplicatedVariants = F)
# optional: filter some samples out of the maf-file
#maf <- filterMaf(maf, tsb = grep("SAMPLE", maf@variants.per.sample$Tumor_Sample_Barcode, value = T, ignore.case = T))
# filter out genes in the maf-file
maf <- filterMaf(maf, genes = "UnknownGene")
maf <- filterMaf(maf, genes = "")
plotmafSummary(maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, fs = 0.8)

oncoplot(maf, showTumorSampleBarcodes = T, SampleNamefontSize = 1, legend_height = 4, 
         legendFontSize = 0.6, barcode_mar = 10, showPct = T, removeNonMutated = F, 
         genes = c("ERG", "STAG2", "BCOR", "CCNB3", "KMTD2", "DUX4", "CIC","EWSR1", "FLI1", "TP53", "FEV"),
       #            "MDN1", "MLLT10", "MUC16", "REG3A", "SLC23A3", "SPINK5", "STMND1",
        #           "USH2A",
        #   "TP53", "STAG2", "MDM2"),
         drawColBar = T, annotationFontSize = 0.6, fontSize = 0.6, 
         titleText = "Mutations in all genes ")


# adding metadata to the Maf object
metadata <- read_xlsx("Metadata.xlsx")
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$TSB3
metadata_maf <- metadata
length(unique(metadata_maf[maf@clinical.data$Tumor_Sample_Barcode, "es-code"]))
maf@clinical.data$Patient_ID <- metadata_maf[maf@clinical.data$Tumor_Sample_Barcode,"M-label"]
maf@clinical.data$FUSION <- metadata_maf[maf@clinical.data$Tumor_Sample_Barcode,"fusion"]
maf@clinical.data$ES_CODE <- metadata_maf[maf@clinical.data$Tumor_Sample_Barcode,"es-code"]
maf@clinical.data$Type <- metadata_maf[maf@clinical.data$Tumor_Sample_Barcode,"biotype"]

Patient_ID_cols <- c("green", "blue", "red", "yellow", "darkblue", 
                "black", "magenta","darkorange", "pink", "gray",
                "orange", "darkgreen", "lightblue", "darkred", "brown",
                "gray30", "gray60", "lightgreen", "coral2", "skyblue", 
                "blue3", "magenta3", "pink3", "green3", "red3", "lightblue3")
FUSION_cols <- c("aquamarine4", "coral3", "darkblue", "yellow", "orange", "magenta", "firebrick")
ES_CODE_cols <- c("green", "blue", "red", "yellow", "darkblue", 
                  "black", "magenta","darkorange", "pink", "gray",
                  "orange", "darkgreen", "lightblue", "darkred", "brown",
                  "aquamarine4", "coral3", "coral1", "seagreen2", "seagreen",
                  "darkgray", "midnightblue", "gray30", "gray20", "coral2",
                  "white", "skyblue", "firebrick", "red3", "blue3", "magenta3", "pink3",
                  "yellow3", "green4", "red4", "lightblue3")
Type_cols <- c("black", "yellow", "coral", "blue")
names(Patient_ID_cols) <-unique(maf@clinical.data$Patient_ID)
names(FUSION_cols) <-unique(maf@clinical.data$FUSION)
names(ES_CODE_cols) <-unique(maf@clinical.data$ES_CODE)
names(Type_cols) <-unique(maf@clinical.data$Type)
maf_filtered_samples <- subsetMaf(maf, tsb = c("M040AAB_tum_wes_ES-055", "M040AAB_org_ES-055", "M660AAA_org_ES-050",
                                               "M660AAA_tum_ES-010", "M660AAA_tum_wes_ES-046", "M660AAA_org_ES-010_2",
                                               "M660AAA_org_ES-046", "M002ZZI_org_ES-002", "M040AAB_org_ES-016",
                                               "M369AAA_org_ES-027", "M616AAA_org_ES-032", "M818AAC_org_ES-057",
                                               "M333AAA_org_ES-071", "M616AAA_org_ES-075", "M040AAB_org_ES-080",
                                               "M003ZZE_pdx_ES-094", "M003ZZF_pdx_ES-096", "M003ZZG_pdx_ES-085",
                                               "M003ZZG_pdx_ES-093", "M002ZZH_org_ES-001", "M057AAA_tum_ES-070",
                                               "M426AAE_tum_ES-102", "M333AAA_tum_ES-042", "M333AAA_tum_wes_ES-071",
                                               "M333AAA_org_ES-042", "M239AAB_tum_ES-041", "M239AAB_org_ES-041",
                                               "M369AAA_tum_ES-027", "M040AAB_tum_ES-016", "M040AAB_tum_wes_ES-080",
                                               "M040AAB_org_ES-016_2", "M818AAC_org_ES-077", "M818AAC_tum_ES-057",
                                               "M818AAC_tum_wes_ES-077", "M081AAD_tum_ES-053", "M081AAD_org_ES-053",
                                               "M616AAA_tum_ES-090", "M616AAA_tum_ES-032", "M616AAA_tum_ES-045", 
                                               "M616AAA_org_ES-045", "M616AAA_org_ES-090", "M003ZZD_org_ES-086",
                                               "M003ZZD_pdx_ES-086"))
oncoplot(maf_filtered_samples, showTumorSampleBarcodes = T, SampleNamefontSize = 1, legend_height = 4, 
         legendFontSize = 0.6, barcode_mar = 6, showPct = T, removeNonMutated = F, 
         genes = c("EWSR1", "FLI1", "ERG", "FEV", "CIC", "DUX4", "BCOR", "KMTD2", "CCNB3", 
                   "STAG2", "TP53", "CREBBP", "MLLT10", "BNC1", "DNMT3B", "ERBB4", "NBN",
                   "TYK2", "MTOR", "TCF7L2", "TNFRSF25", "ARID1A",
                   "CDKN2A", "MDM2", "PTEN", "MYCN", "ALK", "TERT", "IKZF1", "MYC", "CDKN2B", "SERP2", "SETBP1"), 
         drawColBar = T, annotationFontSize = 0.6, fontSize = 0.6, 
         titleText = "Mutations in all genes ", cBioPortal = F, logColBar = T, additionalFeatureCol = "gray99",
         clinicalFeatures = c( "Patient_ID","FUSION", "Type"), bgCol = "gray89", annoBorderCol = "gray90",
         annotationColor = list("Patient_ID"= Patient_ID_cols,
                                "FUSION" = FUSION_cols,
                                "Type" = Type_cols),sortByAnnotation = T,
         sampleOrder = c("M040AAB_tum_ES-016","M040AAB_org_ES-016","M040AAB_org_ES-016_2", 
                         "M040AAB_tum_wes_ES-055", "M040AAB_org_ES-055", "M040AAB_tum_wes_ES-080",
                         "M040AAB_org_ES-080","M660AAA_tum_ES-010", "M660AAA_org_ES-010_2","M660AAA_tum_wes_ES-046", 
                         "M660AAA_org_ES-046", "M660AAA_org_ES-050",
                         "M333AAA_tum_ES-042","M333AAA_org_ES-042","M333AAA_tum_wes_ES-071",
                         "M333AAA_org_ES-071", "M369AAA_tum_ES-027",  "M369AAA_org_ES-027","M239AAB_tum_ES-041",
                         "M239AAB_org_ES-041",     "M002ZZH_org_ES-001",  "M002ZZI_org_ES-002",  
                         "M003ZZG_pdx_ES-085","M003ZZG_pdx_ES-093","M003ZZE_pdx_ES-094",  "M426AAE_tum_ES-102",  
                         "M003ZZF_pdx_ES-096", 
                         "M616AAA_tum_ES-090","M616AAA_org_ES-090",
                         "M616AAA_tum_ES-032", "M616AAA_org_ES-032", "M616AAA_tum_ES-045",
                         "M616AAA_org_ES-045","M616AAA_org_ES-075",
                         "M003ZZD_pdx_ES-086", "M003ZZD_org_ES-086",
                         "M818AAC_tum_ES-057", "M818AAC_org_ES-057","M818AAC_tum_wes_ES-077",
                         "M818AAC_org_ES-077", 
                         "M057AAA_tum_ES-070",
                         "M081AAD_tum_ES-053", 
                         "M081AAD_org_ES-053"),
         tsbToPIDs = data.frame("TSB" = unique(maf_filtered_samples@data$Tumor_Sample_Barcode),
                                "Name" = c("ES-055-T", "ES-055-O", "ES-050-O", "ES-010-T", "ES-046-T", "ES-010-O",
                                           "ES-046-O", "ES-002-O", "ES-016-O", "ES-027-O", "ES-032-O", "ES-057-O",
                                           "ES-071-O", "ES-075-O", "ES-080-O", "ES-094-P", "ES-096-P", "ES-085-P",
                                           "ES-093-P", "ES-001-O", "ES-070-T", "ES-102-T", "ES-042-T", "ES-071-T",
                                           "ES-042-O", "ES-041-T", "ES-041-O", "ES-027-T", "ES-016-T", "ES-080-T",
                                           "ES-016.2-O", "ES-077-O", "ES-057-T", "ES077-T", "ES-053-T", "ES-053-O",
                                           "ES-090-T", "ES-032-T", "ES045-T", "ES-045-O", "ES-090-O", "ES-086-O",
                                           "ES-086-P")))




########### mutational signatures #############
temp <- fread("DNA_analysis/MAF_files_all_sarcomas/all_annotated.maf")
syn <- c("gain", "loss", "gene_fusion", "SNP")
vc <- unique(temp$Variant_Classification)
nonSyn <- setdiff(vc,syn)
maf_temp <- read.maf("DNA_analysis/MAF_files_all_sarcomas/all_annotated.maf", vc_nonSyn = nonSyn)
maf_filtered_with_GL <- subsetMaf(maf_temp,  tsb = c("M040AAB_tum_wes_ES-055", "M040AAB_org_ES-055", "M660AAA_org_ES-050",
                                                                 "M660AAA_tum_ES-010", "M660AAA_tum_wes_ES-046", "M660AAA_org_ES-010_2",
                                                                 "M660AAA_org_ES-046", "M040AAB_org_ES-016",
                                                                 "M369AAA_org_ES-027", "M616AAA_org_ES-032", "M818AAC_org_ES-057",
                                                                 "M333AAA_org_ES-071", "M616AAA_org_ES-075", "M040AAB_org_ES-080",
                                                                 "M057AAA_tum_ES-070","M426AAE_tum_ES-102", "M333AAA_tum_ES-042", 
                                                                 "M333AAA_tum_wes_ES-071", "M333AAA_org_ES-042", "M239AAB_tum_ES-041", 
                                                                 "M239AAB_org_ES-041","M369AAA_tum_ES-027", "M040AAB_tum_ES-016", 
                                                                 "M040AAB_tum_wes_ES-080", "M040AAB_org_ES-016_2", "M818AAC_tum_ES-057",
                                                                 "M818AAC_tum_wes_ES-077", "PMCID974AAI_tum_ES-053", "M081AAD_org_ES-053",
                                                                 "M616AAA_tum_ES-090", "M616AAA_tum_ES-032", "M616AAA_tum_ES-045",
                                                                 "M616AAA_org_ES-045", "M616AAA_org_ES-090"))
maf_filtered_with_GL <- subsetMaf(maf_filtered_with_GL, query = "Variant_Classification %in% nonSyn")

maf_filtered_with_GL.tnm <- trinucleotideMatrix(maf = maf_filtered_with_GL, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", prefix = "chr")
res_test <- estimateSignatures(maf_filtered_with_GL.tnm)
plotCophenetic(res_test)
signatures <- extractSignatures(mat = maf_filtered_with_GL.tnm, n = 8)
test22 <- compareSignatures(signatures, sig_db = "SBS", verbose = TRUE)
plotSignatures(nmfRes = signatures, contributions = T, show_barcodes = T, show_title = T, sig_db = "SBS", absolute = F)#, 
#               patient_order = c())

plotSignatures(nmfRes = signatures, sig_db = "SBS", yaxisLim = NA)
test22$best_match$Signature_1$best_match
pheatmap(signatures$contributions, labels_row = c(test22$best_match$Signature_1$best_match,
                                                                test22$best_match$Signature_2$best_match,
                                                                test22$best_match$Signature_3$best_match,
                                                                test22$best_match$Signature_4$best_match,
                                                                test22$best_match$Signature_5$best_match,
                                                                test22$best_match$Signature_6$best_match,
                                                                test22$best_match$Signature_7$best_match,
                                                                test22$best_match$Signature_8$best_match),
                   labels_col = c("ES-010O", "ES-046O", "ES-050O", "ES-010T", "ES-046T", "ES070T", "ES-042O", "ES-071O",
                                  "ES-042T", "ES-071T", "ES-041O","ES-041T", "ES-027O", "ES-027T", "ES-016O",
                                  "ES-055O", "ES-080O", "ES-016T", "ES-055T", "ES-080T", "ES-057O", "ES-057T", "ES-077T",
                                  "ES-053O", "ES-053T", "ES-032O", "ES-045O", "ES-075O", "ES-090O", "ES-032T", "ES-045T",
                                  "ES-090T"),
                   scale = "column", annotation_row = data.frame("aetiology" = c(test22$best_match$Signature_1$aetiology,
                                                                                  test22$best_match$Signature_2$aetiology,
                                                                                  test22$best_match$Signature_3$aetiology,
                                                                                  test22$best_match$Signature_4$aetiology,
                                                                                  test22$best_match$Signature_5$aetiology,
                                                                                  test22$best_match$Signature_6$aetiology,
                                                                                  test22$best_match$Signature_7$aetiology,
                                                                                  test22$best_match$Signature_8$aetiology),
                                                                 row.names = rownames(signatures$contributions_abs)),
                   color = viridis(100), cluster_rows = T, cluster_cols = T, treeheight_row = 0)




########## Genotype matrix for sample distance tables #############
GT_mat <- genotypeMatrix(maf_filtered_samples, genes = unique(maf_filtered_samples@data$Hugo_Symbol), vafCol = "i_TumorVAF_WU", vafCutoff = c(0,0.6))

GT_mat <- as.data.frame(GT_mat)

GT_mat_num <- GT_mat
for(column in colnames(GT_mat_num)){
  GT_mat_num[,column][which(GT_mat_num[,column] == "None")] <- 0
  GT_mat_num[,column][which(GT_mat_num[,column] == "Het")] <- 1
  GT_mat_num[,column][which(GT_mat_num[,column] == "Hom")] <- 2
  message(paste0(column, "_", sum(as.numeric(GT_mat_num[,column]))))
}

sampleDist <- dist(t(GT_mat_num), method = "binary")
sampleDistMatrix <- as.matrix(sampleDist)
colors <- colorRampPalette(rev(brewer.pal(11, name = "Spectral")))(900)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDist, 
         clustering_distance_cols=sampleDist,
         col=viridis(500, option = "B"),
         main = "sample distance matrix",
         treeheight_col = 0)




