#########################################################################################################################
#
# Gijs van Son
# July 16, 2024
#
# Run pyclone-vi on a sample set
#########################################################################################################################
#!/bin/bash -l

#SBATCH --time=78:0:0
#SBATCH --mem=50G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=g.j.f.vanson@prinsesmaximacentrum.nl

micromamba activate base

pyclone-vi fit -i ${pyclone_input} -o ${library.library_name}.h5 -c 40 -d "beta-binomial -r 100
pyclone-vi write-results-file -i ${library.library_name}.h5 -o ${library.library_name}_pyvi.tsv
