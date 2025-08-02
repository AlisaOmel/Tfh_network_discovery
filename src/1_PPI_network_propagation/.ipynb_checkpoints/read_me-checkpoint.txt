1. example\_hotnet2\_run.slurm --> configuration file for how to run hotnet2
2. get\_network\_prop\_genes.ipynb --> identify significant modules and genes
3. module\_visualization\_fig2A.R --> visualize the significant modules in a random subset of edges.
4. first\_degree\_interactors\_function\_Fig2.ipynb --> identify the first degree interactors from the network and significant genes. Compare across gene sets (like Vinuesa). Make random permutations and gene lists for the upset plots.
5. make\_upset\_plot\_human\_firstdegree\_fig2F.R --> Create the human upset plots seen in Figure 2.
6. forcytoscape\_msigdb\_first\_degree\_network.ipynb --> create the .csvs to input into cytoscape and visualize overlap of the prioritized genes and pathways of interest.
7. steiner\_tree\_extfig2.py --> create shortest path network from a gene list (Vinuesa 2016).
8. run\_steiner\_tree.slurm --> slurm script to run the steiner\_tree\_extfig2.py script
9. steiner\_tree\_extendedgenes\_extfig2.ipynb --> add all interactors of a specific gene to the steiner tree generated network.
