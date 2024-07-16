#########################################################################################################################
#
# Gijs van Son
# July 16, 2024
#
# Construct PyClone-VI input from cnv_facets and vcf-file 
#########################################################################################################################
# load libraries
library(splitstackshape)
library(data.table)

MutectTitan2PyClone <- function(variants, copynum, sample.id, norm.cn = 2) {
  variants <- cbind(variants[,c(1, 2, column.number)], MajorCN = as.numeric(NA), MinorCN = as.numeric(NA))#selects CHROM, POS, AD columns 
  variants[,3] <- unlist(lapply(variants[[3]], function(x){strsplit(x, ":")[[1]][2]}))
  variants <- cSplit(variants, 3, ",") #only takes the first alt_allele on a multiallelic site and ignores the rest
  colnames(variants) <- c("chromosome", "position", "MajorCN", "MinorCN", "ref", "alt")
  copynum2 <- data.frame(CHROM = copynum[,1], START = copynum$POS, END = copynum$INFO, Major_cn = copynum$INFO, Minor_cn = copynum$INFO)
  copynum2$END <- unlist(lapply(copynum2$END, function(x){
    step1 <- strsplit(x, ";")[[1]][3]
    step2 <- strsplit(step1, "ND=")[[1]][2]
    return(step2)}))
  copynum2$Major_cn <- unlist(lapply(copynum2$Major_cn, function(x){
   step1 <- strsplit(x, ";")[[1]]
   step2 <- step1[grep("TCN_EM", step1)]
   step3 <- strsplit(step2, "EM=")[[1]][2]
   return(step3)}))
  copynum2$Minor_cn <-  unlist(lapply(copynum2$Minor_cn, function(x){
    step1 <- strsplit(x, ";")[[1]]
    step2 <- step1[grep("LCN_EM", step1)]
    step3 <- strsplit(step2, "EM=")[[1]][2]
    return(step3)}))
  copynum2$Minor_cn[grep(".", copynum2$Minor_cn, fixed=T)] <- 0
  colnames(copynum2) <- c("CHROM", "START", "END", "major_cn", "minor_cn")
  variants$MajorCN <- 2
  variants$MinorCN <- 1
  for(row in c(1:nrow(variants))){
    if(any(copynum2$CHROM == variants$chromosome[row] & copynum2$START <= variants$position[row] & copynum2$END >= variants$position[row])){
    testrow <- copynum2[copynum2$CHROM == variants$chromosome[row] & copynum2$START <= variants$position[row] & copynum2$END >= variants$position[row],]
    variants[row,3] <- testrow$major_cn[1]
    variants[row,4] <- testrow$minor_cn[1]    
  }}
  id          <- paste(variants$chromosome, variants$position, sep = ':')
  pyclone.tsv <- data.frame(mutation_id = id, sample_id = sample.id, ref_counts = variants$ref, alt_counts = variants$alt,
                            major_cn = variants$MajorCN, minor_cn = variants$MinorCN, normal_cn = norm.cn)
  na.exclude(pyclone.tsv)
}


args <- commandArgs(TRUE)
variants <- fread((args[1]))
copynum <- fread((args[2]))
sample.id <- (args[3])
out.id <- paste0(sample.id, '_pyclone_input.tsv')

input_pyclone <- MutectTitan2PyClone(variants, copynum, sample.id, norm.cn = 2)
write.table(input_pyclone, file=out.id, quote=FALSE, sep='\t', row.names = FALSE)


