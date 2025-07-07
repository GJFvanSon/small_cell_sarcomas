## SUBCLONE ANALYSIS


Scripts to evaluate clusters of mutations that occur over longdituninal samples.

Usage:
# STEP 1

Run CNV facets

change the following parameters in the run_cnv_facets.sh file:
cnv_facets.R \
-t ${PATH to Tumor bam} \
-n ${PATH to Normal bam} \
-vcf ${PATH_TO}00-All.vcf.gz \
-o ${OUTPUT_NAME}

sbatch run_cnv_facets.sh

# STEP 2

convert the cnv facets output to pyclone-VI input using the Rscript

usage:

Rscript make_pyclone_input.R ${VARIANT_VCF_FILE} ${CNV_FACETS_COPYNUMBER_OUT} ${SAMPLE_ID}

This will create an output file called: ${SAMPLE_ID}_pyclone_input.tsv

# STEP3

Run the Run_pyclone.sh script

This will run the following two commands:
pyclone-vi fit -i ${pyclone_input} -o ${library.library_name}.h5 -c 40 -d "beta-binomial -r 100
pyclone-vi write-results-file -i ${library.library_name}.h5 -o ${library.library_name}_pyvi.tsv

This will give you an output tsv file which is processed for figures in R


