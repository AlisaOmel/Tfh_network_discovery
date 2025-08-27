#!/bin/bash
#SBATCH --job-name=RNA_ATAC
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 6-00:00 # Runtime in D-HH:MM
#SBATCH --cpus-per-task=1 # Request that ncpus be allocated per process
#SBATCH --cluster=smp
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=alo90@pitt.edu

module purge
module load taiji/1.2.0
#module load python/ondemand-jupyter-python3.8 gcc/8.2.0
#module load samtools/1.14 sra-toolkit/3.0.2 rgt/0.12.3 taiji/1.2.0

cd /ix/djishnu/Alisa/Tfh/Taiji/ATAC_RNA

export PATH=$PATH:/ix/djishnu/Jane/simar_files/bwa-0.7.17
export PATH=$PATH:/ix/djishnu/Prabal/softwares/samtools-1.16.1

taiji run --config config_atac.yml

echo done
