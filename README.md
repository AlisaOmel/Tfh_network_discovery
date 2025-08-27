# Tfh_network_discovery
Network based approaches for uncovering programs in Tfh differentiation.


# Table of Contents
- [Description](#description)
- [Dependencies](#dependencies)
  - [Packages](#python-packages)
  - [CLI tools](#cli-tools)
- [Citation](#citation)

# Description
This is a set of scripts which perform the following functions:
  1) run a network propagation on a set of gene scores and proccess the output.
  2) sort bulk-RNA data
  3) run a personalized page rank algorithm off of epigenomic and transcriptomic data (Taiji) and process the output.
  4) process scRNA and spatial data and run gene set enrichment analysis of literature curated pathways.
  5) identify pathways of interest from msigDB for both bulk and and scRNA derived data.

These functions are split into folders under Scripts and each contain read_me.txt. 

# Dependencies

Python packages:
- pandas (v 1.2.4)
- numpy
- h5py
- pickle
- random
- itertools
- Seaborn
- scipy
- Networkx (v. 2.5)
- adjustText
- matplotlib
- tqdm
- ast
- csv

R packages:
- Seurat (v5.0.0)
- harmony
- Azimuth
- SingleCellExperiment 
- escape (v0.99.0)
- msigdb (v7.5.1)
- zellkonverter
- BiocParallel
- rogme
- Cairo
- dplyr
- ggpubr
- ggplot2
- qgraph
- viridis
- readr
- ComplexUpset
- pheatmap

CLI tools:
- samtools (v1.16.1)
- bwa (v0.7.17)
- sra-toolkit (v3.0.2) 
- rgt (v0.12.3) 
- taiji (v1.2.0)
- hotnet2 (v1.2.1)

# Citation

