#!/bin/bash

#SBATCH --job-name=steiner_tree
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 6-00:00 # Runtime in D-HH:MM
#SBATCH --cpus-per-task=1 # Request that ncpus be allocated per process
#SBATCH --cluster smp
#SBATCH --mem 60G

python ./steiner_tree.py 
