# Tfh Network Discovery

> Network-based approaches to uncover core programs underlying T follicular helper (Tfh) cell differentiation across bulk, single-cell, and spatial modalities.

---

## Table of Contents
- [Introduction](#introduction)
- [Repository Layout](#repository-layout)
- [Dependencies](#dependencies)
  - [Python](#python)
  - [R](#r)
  - [CLI Tools](#cli-tools)
- [Data & Inputs](#data--inputs)
- [How to Use](#how-to-use)
  - [1) PPI Network Propagation (HotNet2 / Steiner)](#1-ppi-network-propagation-hotnet2--steiner)
  - [2) Bulk RNA Processing](#2-bulk-rna-processing)
  - [3) Taiji Network Propagation](#3-taiji-network-propagation)
  - [4) scRNA & Spatial Processing + GSEA](#4-scrna--spatial-processing--gsea)
  - [5) Pathway Discovery (Bulk & scRNA)](#5-pathway-discovery-bulk--scrna)
- [Linking Inputs → Code → Outputs](#linking-inputs--code--outputs)
- [Reproducibility Notes](#reproducibility-notes)
- [Citation](#citation)

---

## Introduction
This repository hosts a modular analysis pipeline that (i) propagates transcriptomic/epigenomic signals through protein–protein interaction (PPI) and gene‑regulatory networks, (ii) integrates bulk RNA‑seq, scRNA‑seq, and spatial transcriptomics, and (iii) performs pathway discovery using MSigDB gene sets. Outputs include prioritized genes/modules, enrichment summaries, and publication‑ready tables/figures.

Each stage is self‑contained under `Scripts/` with minimal coupling, enabling you to run individual steps or the full workflow.

---

## Repository Layout
```
.
├── Input_data/                 # Primary inputs (PPI graph, labels, Taiji outputs, curated sets)
├── Sample_outputs/             # Small example outputs to validate runs
├── Scripts/                    # Analysis code organized by stage
│   ├── 1_PPI_network_propagation/
│   ├── 2_bulk_RNA_processing/
│   ├── 3_Taiji_network_propagation/
│   ├── 4_scRNA_processing_and_analysis/
│   ├── 5_pathway_discovery_bulk/
│   └── 5_pathway_discovery_scRNA/
└── supplement_tables/          # Full-size result tables used in the paper/SI
```

---

## Dependencies

### Python
- pandas (≥ 1.2.4)
- numpy (≥ 1.24.3)
- h5py (≥ 3.9.0)
- scipy (≥ 1.11.1)
- networkx (≥ 2.5)
- matplotlib (≥ 3.7.2)
- seaborn (≥ 0.12.2)
- adjustText (≥ 1.3.0)
- tqdm (≥ 4.65.0)

### R
- Seurat (≥ 5.3.0), harmony (≥ 1.2.3), Azimuth (≥ 0.5.0)
- SingleCellExperiment(≥ 1.30.1), zellkonverter(≥ 1.18.0), BiocParallel (≥ 1.42.1)
- escape (≥ 0.99.0)
- msigdb (≥ 7.5.1)
- rogme (≥ 0.2.1), Cairo (≥ 1.6.2), dplyr (≥ 1.1.4), ggpubr (≥ 0.6.1), ggplot2(≥ 3.5.2), qgraph (≥ 1.9.8), viridis (≥ 0.6.5), readr(≥ 2.1.5)
- ComplexUpset(≥ 1.3.3), pheatmap (≥ 1.0.13)

### CLI Tools
- samtools (≥ 1.16.1)
- bwa (≥ 0.7.17)
- SRA Toolkit (≥ 3.0.2)
- RGT (≥ 0.12.3)
- Taiji (≥ 1.2.0)
- HotNet2 (≥ 1.2.1)

---

## Data & Inputs
Key inputs live under `Input_data/`:
- **PPI graph / propagation inputs**: `HomoSapiens_binary_co_complex_Feb2023_1_ppr_0.4.h5`, curated `.sif` edges, and node attributes.
- **Curated gene sets**: `literature_sets/` and various `*_unique_genes*.pkl/csv`.
- **Bulk & Taiji processed data**: `processed_Taiji/*.csv`, matrices and ranked files.
- **Labels for scRNA comparisons**: `labels_GUT_*_combo.rds`.

Large binary assets are stored via Git LFS.

---

## How to Use

### 1) PPI Network Propagation (HotNet2 / Steiner)
**Where**: `Scripts/1_PPI_network_propagation/`

- `HotNet2.py` & `example_hotnet2_run.slurm`: run diffusion/propagation over the PPI graph.
- `steiner_tree.py` / `steiner_tree_extfig2.py` & `run_steiner_tree.slurm`: compute Steiner sub‑networks connecting seeds.
- Notebooks (`get_network_prop_genes.ipynb`, `steiner_tree_extendedgenes_extfig2.ipynb`) for exploration/exports.
- **Inputs**: PPI HDF5 (`Input_data/HomoSapiens_*.h5`), seed/score lists (e.g., `Input_data/pps_unique_genes.pkl`), `.sif` files (`extended_vinuesa_*.sif`).
- **Example outputs**: `Sample_outputs/pps_protein_pairs_sig_genes_df.csv`; SI tables under `supplement_tables/*overlaps_df*.csv`.

**Run examples**
```bash
# Slurm example (HotNet2)
sbatch Scripts/1_PPI_network_propagation/example_hotnet2_run.slurm

# Slurm example (Steiner)
sbatch Scripts/1_PPI_network_propagation/run_steiner_tree.slurm
```

### 2) Bulk RNA Processing
**Where**: `Scripts/2_bulk_RNA_processing/`

- `sort_bulk_RNA_forTaiji.ipynb`, `logFCrna_calculation.ipynb`.
- **Inputs**: `Input_data/Run_313.*normalized_data_matrix.tsv` and related bulk RNA matrices.
- **Outputs**: separated `.tsv` per state under `Sample_outputs/sorted_bulk/separated_tsv/` and logFC summaries.

### 3) Taiji Network Propagation
**Where**: `Scripts/3_Taiji_network_propagation/`

- Configs: `config_atac.yml`, `mod_input_atac.yml`, plus `config_atac.slurm`.
- `processing_Taiji_output.ipynb` post‑processes Taiji results.
- **Inputs**: ranked/processed CSVs in `Input_data/processed_Taiji/`.
- **Outputs**: prioritized Taiji gene lists (e.g., `Sample_outputs/Taiji_nodot_genes.pkl`) and combined tables in `supplement_tables/`.

**Run example**
```bash
sbatch Scripts/3_Taiji_network_propagation/config_atac.slurm
```

### 4) scRNA & Spatial Processing + GSEA
**Where**: `Scripts/4_scRNA_processing_and_analysis/`

- Human/mouse processing scripts (`human_scRNA_gut_processing.R`, `mouse_scRNA_LN_analysis_4CDE.R`), spatial analyses (`human_spatial_scRNA_analysis_3H.R`).
- `example_ssGSEA_fromGeneSets_plotting_fig3DE.R` to compute ssGSEA and generate figures.
- **Inputs**: label RDS files under `Input_data/labels_*.rds`; curated gene sets.
- **Outputs**: enrichment summaries (e.g., `supplement_tables/hallmark_ssGSEA_*summary*.csv`, `pid_ssGSEA_*summary*.csv`).

### 5) Pathway Discovery (Bulk & scRNA)
**Where**: `Scripts/5_pathway_discovery_bulk/` and `Scripts/5_pathway_discovery_scRNA/`

- Bulk: `example_permutations_for_bulk_msigdb_analysis.R`, `make_volcano_plot.R`, Slurm wrapper `example_run_permutations_for_volcanoplot_fig5.slurm`.
- scRNA: `example_scGSEA_msigdb_fig5.R`, `pathway_boxplots_5CFI.R`.
- **Inputs**: MSigDB catalogs, seed sets, and modality‑specific score tables.
- **Outputs**: volcano tables (e.g., `supplement_tables/patternb_taiji_{HALLMARK|KEGG|PID}_1000_volc_plot_thresh_qval.csv`) and figures.

---

## Linking Inputs → Code → Outputs

| Stage | Primary Scripts | Key Inputs (examples) | Expected Outputs (examples) |
|---|---|---|---|
| 1. PPI propagation | `HotNet2.py`, `steiner_tree.py`, `run_steiner_tree.slurm` | `Input_data/HomoSapiens_*0.4.h5`, `Input_data/pps_unique_genes.pkl`, `.sif` and node attributes | `Sample_outputs/pps_protein_pairs_sig_genes_df.csv`, `supplement_tables/random_*overlaps_df*.csv` |
| 2. Bulk RNA | `logFCrna_calculation.ipynb`, `sort_bulk_RNA_forTaiji.ipynb` | `Input_data/Run_313.*normalized_data_matrix.tsv` | `Sample_outputs/sorted_bulk/separated_tsv/*.tsv`, logFC summaries |
| 3. Taiji | `config_atac.yml`, `processing_Taiji_output.ipynb` | `Input_data/processed_Taiji/*.csv` | `Sample_outputs/Taiji_nodot_genes.pkl`, `supplement_tables/final_genes_all_sets_ALLGENES.csv` |
| 4. scRNA & spatial + GSEA | `human_scRNA_gut_processing.R`, `mouse_scRNA_LN_analysis_4CDE.R`, `example_ssGSEA_fromGeneSets_plotting_fig3DE.R` | `Input_data/labels_GUT_*_combo.rds`, curated sets in `Input_data/literature_sets/` | `supplement_tables/*ssGSEA*summary*.csv`, `supplement_tables/cliffs_delta_summary_*.csv` |
| 5. Pathway discovery | `example_permutations_for_bulk_msigdb_analysis.R`, `make_volcano_plot.R`, `example_scGSEA_msigdb_fig5.R` | MSigDB catalogs; outputs from stages 1–4 | `supplement_tables/patternb_taiji_*_volc_plot_thresh_qval.csv`, boxplots/figures |

---

## Reproducibility Notes
- Large inputs are tracked with **Git LFS**. Ensure `git lfs install` before cloning.
- Some scripts assume **Slurm** availability; adapt to your scheduler or run locally where feasible.
- Randomization/permutation steps set seeds inside scripts; for exact reproduction, also control global RNG (R: `set.seed()`, Python: `numpy.random.seed()`).

---

## Citation
If you use this code or its outputs, please cite the associated manuscript and tools/databases referenced in the scripts (Taiji, HotNet2, MSigDB, Seurat, etc.). Add the full citation here once available.

