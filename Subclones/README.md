# SUBCLONE ANALYSIS


Scripts to evaluate clusters of mutations that occur over longdituninal samples.

Usage:

Run CNV facets

change the following parameters in the run_cnv_facets.sh file:
cnv_facets.R \
-t ${PATH to Tumor bam} \
-n ${PATH to Normal bam} \
-vcf ${PATH_TO}00-All.vcf.gz \
-o ${OUTPUT_NAME}

sbatch run_cnv_facets.sh

