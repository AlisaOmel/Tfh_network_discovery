#!/bin/bash

#SBATCH --job-name=hallmark_1000perm
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 2-00:00 # Runtime in D-HH:MM
#SBATCH --cpus-per-task=1 # Request that ncpus be allocated per process
#SBATCH --cluster htc
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=alo90@pitt.edu
#SBATCH --mem 60G
#SBATCH --output '/ix3/djishnu/Alisa/Tfh/ForPaper/Scripts/PPI_network_propagation/hallmark_1000perm.out'

cd /ix3/djishnu/Alisa/Tfh/ForPaper/Scripts/PPI_network_propagation
  
Rscript example_permutations_for_bulk_msigdb_analysis.R

echo done