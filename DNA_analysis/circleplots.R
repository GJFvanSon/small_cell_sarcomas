############################################################################################
#            Creating circleplots with copynumber, fusion and mutation info                #
#                                                                                          #
#                                   Gijs van Son                                           #          
#                                   June 11 2024                                           #
#                         g.j.f.vanson@prinsesmaximacentrum.nl                             #                            
#                                                                                          #
############################################################################################
########### setting the environment #############
library(maftools)
library(data.table)
library(circlize)
library(NMF)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)

setwd("/Users/gijsvanson/OneDrive - Prinses Maxima Centrum/Femke/CIC_DUX4/")

########## defining the circleplot function #############
make_circos_plot <- function(snv_in, ratiofile_in, sv_in = NULL, genome = "hg38",
                             samplename, outputfile, WES = F, SEX="M", SV_type = "Manta",
                             jitter=F, coding_snps_only=F, color = "orange", CNV_colors = c("skyblue2", "gray77", "firebrick2")){
  # get all the SNVs from the file and make them into a table
  #check variables
  if(is.null(snv_in)){
    return("ERROR: No SNP-file found, please specify snv_in = yourfile.vcf")}
  if(is.null(ratiofile_in)){
    return("ERROR: No CNV-file found, please specify ratiofile_in = yourfile.csv")}
  if(genome != "hg38"){return("ERROR: Only hg38 is supported so far!")}
  if(is.null(samplename)){return("ERROR:Please specify a sample name by samplename = yourname")}
  if(isFALSE(SV_type %in% c("Manta", "Fusion"))){
    return("ERROR: Choose between Manta or Fusion for SV_type")}
  if(isFALSE(SEX %in% c("M", "F"))){return("ERROR: sex should be either M(ale) or F(emale)")}
  
  circos.clear()
  snv_file_info <- file.info(snv_in)
  test1 <- system(paste0("cat ", snv_in, " | wc -l"), intern = T)
  test2 <- system(paste0("cat ", snv_in, " | grep ^# | wc -l"), intern = T)
  if (test1 != test2) {
    snv_check <- "full"
    snv_tmp <- read.table(snv_in,comment.char="#")
    if(coding_snps_only){
      if(length(unique(grep("MODERATE", snv_tmp$V8), grep("HIGH", snv_tmp$V8)))==0){
        message("No SNPs found in coding regions")
        snv_check <- "empty"
      }else{
        snv_tmp <- snv_tmp[unique(grep("MODERATE", snv_tmp$V8), grep("HIGH", snv_tmp$V8)),]
      }}
    snv <- as.data.frame(cbind(paste0("chr",as.character(snv_tmp[,1])), snv_tmp[,2], snv_tmp[,2]+1, unlist(lapply(snv_tmp$V10, function(x){strsplit(x,":")[[1]][12]}))))
    colnames(snv) <- c("CHROM", "START", "END", "VAF")
    if(any(snv$CHROM %in% c("chrMT", grep("chrUn", snv$CHROM, value = T)))){
      snv <- snv[-which(snv$CHROM %in% c("chrMT", grep("chrUn", snv$CHROM, value = T))),]
    }
    snv$CHROM <- as.character(snv$CHROM)
    snv$START <- as.numeric(snv$START)
    snv$END <- as.numeric(snv$END)
    snv$VAF <- as.numeric(snv$VAF)
    if(snv_check == "empty"){snv = NULL}
  }
  if (test1 == test2) {
    snv <- NULL
    message("No SNPs found in chromosomes 1 to 22 + X and Y regions")
  }
  
  # get the copynumbers from the ratiofiles
  cnv <- read.table(ratiofile_in, header = T)
  if(WES){
    cnv2 <- data.frame("Chromosome" = cnv$Chromosome, "Start" = cnv$Start, "End" = unlist(lapply(cnv$Gene, function(x){strsplit(x, "-")[[1]][2]})), 
                       "Ratio" = cnv$Ratio, "CopyNumber" = cnv$CopyNumber)
    # correct for sex differences in chromosome copy numbers
    if (SEX == "F") {
      cnv2$CopyNumber[which(cnv2$Chromosome == "Y")] <- 2
    }
    if (SEX == "M") {
      cnv2$CopyNumber[which(cnv2$Chromosome == "X")] <- cnv2$CopyNumber[which(cnv2$Chromosome == "X")]*2
      cnv2$CopyNumber[which(cnv2$Chromosome == "Y")] <- cnv2$CopyNumber[which(cnv2$Chromosome == "Y")]*2
    }
  }
  if (!WES) {
    cnv2 <- data.frame("Chromosome" = cnv$Chromosome, "Start" = cnv$Start, "End" = unlist(lapply(cnv$Start, function(x){x+999})), 
                       "Ratio" = cnv$Ratio, "CopyNumber" = cnv$CopyNumber)
    test <- lapply(seq(1,dim(cnv2)[[1]], 50)[-length(seq(1,dim(cnv2)[[1]], 50))], function(x){
      if(all(cnv2$Chromosome[c(x,x+49)] == cnv$Chromosome[x])){
        chr <- cnv2$Chromosome[x]
        start <- cnv2$Start[x]
        end <- cnv2$End[x+49]
        ratio <- median(cnv2$Ratio[c(x:x+49)])
        # correct for sex differences in chromosome copy numbers
        if (chr == "X") {
          if (SEX == "M") {
            copynumber <- median(cnv2$CopyNumber[c(x:x+49)])*2
          }
          if (SEX == "F"){
            copynumber <- median(cnv2$CopyNumber[c(x:x+49)])
          }
        }
        if (chr == "Y"){
          if(SEX == "M"){
            copynumber <- median(cnv2$CopyNumber[c(x:x+49)])*2
          }
          if(SEX == "F"){
            copynumber <- 2
          }}
        if (chr %in% c(1:22)) {
          copynumber <- median(cnv2$CopyNumber[c(x:x+49)])
        }
        tmp <- data.frame("Chromosome" = chr, "Start" = start, "End" = end, 
                          "Ratio" = ratio, "CopyNumber" = copynumber)
        return(tmp)
      }
      
    })
    test2 <- test[-which(unlist(lapply(test, function(x){length(x)})) != 5)]
    test2 <- do.call("rbind",test2)
    cnv2 <- test2
  }
  cnv3 <- cnv2
  cnv3[,1]=paste("chr",as.character(cnv3[,1]),sep="")
  cnv3$CopyNumber[which(cnv3$CopyNumber > 4)] <- 4
  # make the structual variants table ready for use
  if(SV_type == "Fusion"){
    if (!is.null(sv_in)) {
      sv <- fread(sv_in)
      sv <- as.data.frame(sv)
      all_chrom <- unique(unlist(lapply(c(1:length(sv$LeftGene)), function(x){
        return(c(strsplit( sv$LeftBreakpoint[x], ":")[[1]][1], strsplit(sv$RightBreakpoint[x], ":")[[1]][1]))  
      })))
      if(length(all_chrom) > 0){
        svTable <- data.frame(unlist(lapply(c(1:length(sv$LeftGene)), function(x){
          strsplit(sv$LeftBreakpoint[x], ":")[[1]][1]})),
          unlist(lapply(c(1:length(sv$LeftGene)), function(x){
            strsplit(sv$LeftBreakpoint[x], ":")[[1]][2]
          })),
          unlist(lapply(c(1:length(sv$LeftGene)), function(x){
            strsplit(sv$RightBreakpoint[x], ":")[[1]][1]})),
          unlist(lapply(c(1:length(sv$LeftGene)), function(x){
            strsplit(sv$RightBreakpoint[x], ":")[[1]][2]
          })))
        
      }
      # Filter out all variants that map partly or completely to unknown chromosomes
      if(any(svTable[,1] %in% grep("chrUn", all_chrom, value = T))){
        if(!all(svTable[,1] %in% grep("chrUn", all_chrom, value = T))){
          if(length(grep("Un", svTable[,1])) > 0){
            svTable <- svTable[-grep("Un", svTable[,1]),]}
        }}
      if(any(svTable[,1] == "chrMT")){
        if(!all(svTable[,1] == "chrMT")){
          if(length(grep("MT", svTable[,1])) > 0){
            svTable <- svTable[-grep("MT", svTable[,1]),]}
        }}
      if (any(svTable[,1] %in% c("chrMT", grep("chrUn", all_chrom, value = T)))) {
        sv_in <- NULL
      }
      if(any(svTable[,3] %in% grep("chrUn", all_chrom, value = T))){
        if(!all(svTable[,3] %in% grep("chrUn", all_chrom, value = T))){
          if(length(grep("Un", svTable[,3])) > 0){
            svTable <- svTable[-grep("Un", svTable[,3]),]}
        }}
      if(any(svTable[,3] == "chrMT")){
        if(!all(svTable[,3] == "chrMT")){
          if(length(grep("MT", svTable[,3])) > 0){
            svTable <- svTable[-grep("MT", svTable[,3]),]}
        }}
      if (any(svTable[,3] %in% c("chrMT", grep("chrUn", all_chrom, value = T)))) {
        sv_in <- NULL
      }}}
  if(SV_type == "Manta"){
    if (!is.null(sv_in)) {
      sv <- read.table(sv_in,comment.char="#")
      colnames(sv) <- c("CHROM", "START", "TYPE", "REF", "ALT", "V6", "FILTER", "INFO", "GT_legend", "GT")
      sv <- sv[sv$FILTER == "PASS",]
      sv <- sv[which(sv$V6 > 998),]
      sv <- sv[which(unlist(lapply(sv$GT, function(x){
        temp <- strsplit(x, ":")[[1]][6]
        temp2 <- strsplit(temp, ",")[[1]][2]
        return(temp2)})) > 20),]
      all_chrom <- unique(unlist(c(lapply(sv$CHROM, function(x){paste0("chr",x)}), "chrUnxx")))
      if(any(unlist(lapply(sv$TYPE, function(x){strsplit(x, ":")[[1]][1]})) == "MantaBND")){
        sv_bnd <- sv[grep("BND", sv$TYPE),]
        svTable=data.frame(paste0("chr", sv_bnd[,1]), as.numeric(sv_bnd[,2]), 
                           paste0("chr",as.character(unlist(lapply(sv_bnd$ALT, function(x){strsplit(strsplit(x,"]")[[1]][2],":")[[1]][1]})))),
                           as.numeric(unlist(lapply(sv_bnd$ALT, function(x){strsplit(strsplit(x,"]")[[1]][2],":")[[1]][2]}))))
        svTable <- svTable[-which(is.na(svTable[,4])),]
      }
      # Filter out all variants that map partly or completely to unknown chromosomes
      if(any(svTable[,1] %in% grep("chrUn", all_chrom, value = T))){
        if(!all(svTable[,1] %in% grep("chrUn", all_chrom, value = T))){
          if(length(grep("Un", svTable[,1])) > 0){
            svTable <- svTable[-grep("Un", svTable[,1]),]}
        }}
      if(any(svTable[,1] == "chrMT")){
        if(!all(svTable[,1] == "chrMT")){
          if(length(grep("MT", svTable[,1])) > 0){
            svTable <- svTable[-grep("MT", svTable[,1]),]}
        }}
      if (any(svTable[,1] %in% c("chrMT", grep("chrUn", all_chrom, value = T)))) {
        sv_in <- NULL
      }
      if(any(svTable[,3] %in% grep("chrUn", all_chrom, value = T))){
        if(!all(svTable[,3] %in% grep("chrUn", all_chrom, value = T))){
          if(length(grep("Un", svTable[,3])) > 0){
            svTable <- svTable[-grep("Un", svTable[,3]),]}
        }}
      if(any(svTable[,3] == "chrMT")){
        if(!all(svTable[,3] == "chrMT")){
          if(length(grep("MT", svTable[,3])) > 0){
            svTable <- svTable[-grep("MT", svTable[,3]),]}
        }}
      if (any(svTable[,3] %in% c("chrMT", grep("chrUn", all_chrom, value = T)))) {
        sv_in <- NULL
      }}}
  # make the actual plots
  circos.clear()
  circos.par("start.degree" = 90, "track.height" = 0.13, cell.padding = c(0, 0, 0, 0) )
  if(SEX == "M"){
    circos.initializeWithIdeogram(species = "hg38")}
  if(SEX == "F"){
    circos.initializeWithIdeogram(species = "hg38", chromosome.index = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                                                                         "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
                                                                         "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
                                                                         "chrX"))
  }
  
  col_fun = colorRamp2(breaks = c(1.5,2,2.5), colors = CNV_colors)
  title(samplename)
  if (isFALSE(is.null(snv))) {
    circos.genomicTrackPlotRegion(snv,stack = F ,ylim = c(0,1),panel.fun = function(region, value, ...) {
      circos.genomicPoints(region, track.index = 3,value ,cex = 0.5, pch = 16,col=color , ...)
    })
  }
  if (is.null(snv)) {
    snv <- data.frame("col1" = NA, "col2" = NA, "col3" = NA, "col4" = NA)
    circos.genomicTrackPlotRegion(snv,stack = F,ylim = c(0,1), panel.fun = function(region, value, ...) {
      circos.genomicPoints(region, value, track.index = 3,cex = 0.05, pch = 9,col=color , ...)
    })
  }
  if (jitter){ 
    circos.genomicTrackPlotRegion(cnv3[,c(1,2,3,5)], stack = F, ylim = c(0, 4.5),panel.fun = function(region, value, ...){
      circos.genomicRect(region, value, track.index = 4, col = col_fun(value[[1]]), border = NA)
      circos.genomicLines(region, value,lwd = 1,track.index = 4,baseline = 0, ...)
    })}else{
      circos.genomicTrackPlotRegion(cnv3[,c(1,2,3,5)], stack = F, ylim = c(0, 4.5),panel.fun = function(region, value, ...){
        circos.genomicRect(region, value, track.index = 4, col = col_fun(value[[1]]), border = NA)
      })}
  
  if (is.null(sv_in)){message("no structual variants found")}
  if (!is.null(sv_in)) {
    if(!is.null(svTable)){
      # svTable <- svTable[-which(svTable$paste0..chr...sv_bnd...1.. == svTable$paste0..chr...as.character.unlist.lapply.sv_bnd.ALT..function.x...),]
      for (i in 1:dim(svTable)[1]) {
        circos.link(as.character(svTable[i,1]), as.numeric(svTable[i,2]), as.character(svTable[i,3]), as.numeric(svTable[i,4]))
      }
      legend(0.63,0.845,legend="SV-TRANSLOCATION",col="black",lty=1,cex=0.75,lwd=1.2,bty='n')
    }}
  if(jitter){
    if(coding_snps_only){
      legend(0.7, 1.02, legend = c("Coding SNVs", "CNV-DUPLICATION","CNV-DELETION"), col = c(color, CNV_colors[c(1,3)]),pch=c(16,15,15),cex=0.75,title="Tracks:",bty='n')
    }else{
      legend(0.7,1.02,legend=c("SNV", "CNV-DUPLICATION","CNV-DELETION"),col=c(color,CNV_colors[c(1,3)]),pch=c(16,15,15),cex=0.75,title="Tracks:",bty='n')}
  }
  if(!is.null(outputfile)){
    dev.copy2pdf(file = outputfile)
  }
  #clear up mess
  circos.clear()
}

################ calling the function to create circleplots ######################
make_circos_plot(snv_in = "DNA_analysis/SMURF_output/ES-016.SMuRF.filtered.vcf",
                 sv_in = "RNA_fusion_outputs/PMOBM000ABK_PMCRZ057HFK_RNA-Seq.star-fusion_predicted.annotated.filtered.tsv", SEX = "M",
                 ratiofile_in = "DNA_analysis/CNV_files/ES-016_dedup.bam_ratio.txt", outputfile = "Figures_output/circle_es-016_org",
                 samplename = "ES-016_org", SV_type = "Fusion", jitter=T, color = "midnightblue")


