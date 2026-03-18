# Cohesin loop extrusion analysis

This repository contains the scripts used to analyze and generate the major figures for the manuscript:

**“Cohesin extrudes DNA unidirectionally through two modes of action in human cells.”**

## Repository structure

The repository is organized by analysis stage:

- **01_preprocessing/**  
  Scripts for defining genomic features, including CTCF motifs, convergent loops, and cohesin loading sites.

- **02_basic_analysis/**  
  Core statistical and descriptive analyses, including enrichment, correlation, and quantitative characterization of chromatin features.

- **03_pileup_analysis/**  
  Aggregated analyses centered on loop anchors and cohesin loading sites, including normalization and pileup of contact maps.

- **04_simulation/**  
  Scripts for downstream analysis of simulated contact maps.  
  The simulation framework is available at: https://doi.org/10.5281/zenodo.19093829

- **05_visualization/**  
  Scripts for generating figures, including contact map visualization and multi-track genomic plots.

## Contributors

Li Yang  
Linghan Jiang  
Xiaotao Wang  

## Citation

If you use this code, please cite:

Wang, P., Meng, L., Jiang, L., Yang, L., Huang, J., Yu, T., Chai, H., Kim, M., Wang, X., Ruan, Y. 
*Cohesin extrudes DNA unidirectionally through two modes of action in human cells.* 2026