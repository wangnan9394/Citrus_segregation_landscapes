# Code Repository Description
This repository contains all analytical workflows and scripts developed for the study of segregation distortion uncovered during the evolutionary transition from apomixis to sexual reproduction in *Citrus hindsii*. The code supports comprehensive genomic and genetic analyses across three core research dimensions, enabling reproducible reconstruction of haplotype dynamics, recombination landscapes, segregation distortion modeling, and pangenome architecture.

## Core Analytical Workflows
The repository is organized into three key modules, corresponding to the main research directions of the project:

### 1. Genome Assembly and Pangenome Construction
Workflows for de novo genome assembly and haplotype-resolved pangenome construction, including pipelines for high-quality phased genome assembly, chromatin interaction (Hi-C) scaffolding, and graph-based pangenome generation using both minigraph-cactus and PGGB frameworks. Downstream analyses such as structural variation detection and genome annotation are also included.

### 2. Haplotype-Specific Marker Identification and Progeny Haplotype Shuffling Map Construction
Integrated pipelines for identifying haplotype-specific markers and reconstructing the genome-wide haplotype shuffling landscape of hybrid progeny. Key steps include diagnostic marker calling, high-resolution recombination breakpoint detection, bin partitioning based on recombination events, and haplotype origin inference across the full genetic map.

### 3. Exponential Decay Model Construction for Segregation Distortion
Complete implementation of the exponential decay model for quantitative dissection of segregation distortion, covering three interconnected analytical sub-modules:
- Single-locus model fitting and parameter estimation
- Multi-locus superposition modeling to simulate combined effects of independent causal variants
- Model calibration, fitting, and causal locus mapping based on empirical population genomic data

## Notes
All scripts are modular, well-documented, and designed for reproducible large-scale genomic analysis, supporting the core conclusions of haplotype-specific segregation distortion, recombination dynamics, and the genomic basis underlying the reproductive mode transition in *Citrus hindsii*.
