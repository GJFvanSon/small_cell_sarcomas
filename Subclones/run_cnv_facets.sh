#########################################################################################################################
#
# Gijs van Son
# July 16, 2024
#
# Run CNV_facets to obtain minor and major copy numbers from genomic data 
#########################################################################################################################


#!/bin/bash -l

#SBATCH --time=78:0:0
#SBATCH --mem=20G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=g.j.f.vanson@prinsesmaximacentrum.nl

micromamba activate cnv_facets

cnv_facets.R \
-t /home/clonevo-gvson/Projects/Ewing_sarcoma/BAMS/PMLBM000FOA-X_dedup.bam \
-n /home/clonevo-gvson/Projects/Ewing_sarcoma/BAMS/PMLBM000FMX-X_dedup.bam \
-vcf /home/clonevo-gvson/General/00-All.vcf.gz \
-o /home/clonevo-gvson/Projects/Ewing_sarcoma/ES-080T


