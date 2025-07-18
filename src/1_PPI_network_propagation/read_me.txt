example_hotnet2_run.slurm --> configuration file for how to run hotnet2
get_network_prop_genes.ipynb --> identify significant modules and genes
module_visualization_fig2A.R --> visualize the significant modules in a random subset of edges.
first_degree_interactors_function_Fig2.ipynb --> identify the first degree interactors from the network and significant genes. Compare across gene sets (like Vinuesa). Make random permutations and gene lists for the upset plots.
make_upset_plot_human_firstdegree_fig2F.R --> Create the human upset plots seen in Figure 2. 
forcytoscape_msigdb_first_degree_network.ipynb --> create the .csvs to input into cytoscape and visualize overlap of the prioritized genes and pathways of interest.
steiner_tree_extfig2.py --> create shortest path network from a gene list (Vinuesa 2016).
run_steiner_tree.slurm --> slurm script to run the steiner_tree_extfig2.py script
steiner_tree_extendedgenes_extfig2.ipynb --> add all interactors of a specific gene to the steiner tree generated network.

 

