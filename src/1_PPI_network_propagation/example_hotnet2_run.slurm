#!/bin/bash

#SBATCH -t 6-00:00 # Runtime in D-HH:MM
#SBATCH --job-name=B_hotnet_TFH
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=alo90@pitt.edu
#SBATCH --error=/ix/djishnu/Alisa/Tfh/Network_analysis/results/patternB/hotnet.err

module purge
module load python/anaconda2.7-4.2.0

# Set location of HotNet2, number of cores, number of network permutations, and number of heat permutations.
#hotnet2=..
num_cores=-1

python /ix/djishnu/TFH/hotnet2-master/HotNet2.py \
-nf  /ix/djishnu/Priyamvada/Immport/Network_analysis/Networks/Networks_to_run_hotnet/HomoSapiens_binary_co_complex_Feb2023_1_ppr_0.4.h5 \
-pnp /ix/djishnu/Priyamvada/Immport/Network_analysis/Networks/Networks_to_run_hotnet/permuted_networks/HomoSapiens_binary_co_complex_Feb2023_1_ppr_0.4_##NUM##.h5 \
-hf  /ix/djishnu/Alisa/Tfh/Network_analysis/data/Pattern_1234_Genes_Scores.txt\
-np 500 \
-hp 500 \
-o /ix/djishnu/Alisa/Tfh/Network_analysis/results/patternB \
-c  $num_cores
