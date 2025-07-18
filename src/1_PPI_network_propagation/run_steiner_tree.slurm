#!/bin/bash

#SBATCH --job-name=steiner_tree
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 6-00:00 # Runtime in D-HH:MM
#SBATCH --cpus-per-task=1 # Request that ncpus be allocated per process
#SBATCH --cluster smp
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=alo90@pitt.edu
#SBATCH --mem 60G

python /ix/djishnu/Alisa/Tfh/ForPaper/steiner_tree.py 
