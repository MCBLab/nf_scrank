# NF_scRank 🧬

**MCBLab/NF_scRank** is a scalable Nextflow pipeline designed to infer Gene Regulatory Networks (GRNs) and calculate single-cell expression ranking perturbation scores using the GENIE3 algorithm. 

Previous GRN tools are often difficult to scale for large Single-Cell RNA-seq (scRNA-seq) datasets. In this context, `NF_scRank` was built to enable high-throughput perturbation scoring in a user-friendly, parallelized, and computationally effective way. The pipeline uses Singularity containers, making installation trivial and results highly reproducible across high-performance computing (HPC) environments.

## Pipeline Summary

The workflow executes the following core modules:

### 1. Object Parsing and Downsampling (`DOWNSAMPLE`)
This is the initial step of the process. It ingests a fully processed Seurat object (`.rds`) and identifies the user-defined metadata column containing the cell identities (e.g., cell types or clones). To ensure statistical robustness and equitable GRN inference, it randomly downsamples the cells from each identity to a specified maximum number (`--n_cells`), balancing the computational load.

### 2. Expression Matrix Extraction (`SPLIT_MATRICES`)
The pipeline isolates the RNA assay from the downsampled Seurat object and generates independent `.csv` expression matrices for each cellular identity. This step isolates the transcriptional state of each group.

### 3. GENIE3 Network Inference (`GENIE3_NETWORK`)
This is the heavy-lifting computational core. For each isolated expression matrix, the pipeline runs [GENIE3](https://bioconductor.org/packages/release/bioc/html/GENIE3.html) (GEne Network Inference with Ensemble of trees). It infers a co-expression regulatory network where nodes are genes and edges represent the "weight" of the regulatory link between them, returning a comprehensive matrix of regulatory interactions per cell type.

### 4. Perturbation Scoring (`SCRANK_SCORE`)
Using the list of target genes (`--target`) provided by the user, this module extracts the specific regulatory weight of the targets from the GENIE3 output. It calculates the perturbation score, which reflects how much the network relies on the specific target gene within that specific cell state.

### 5. Consolidate Results (`MERGE`)
This final step collects the perturbation scores from all parallel GENIE3 tasks and merges them into a single, clean text file, ready for downstream visualization.

## Quick Start
1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html) (`>=22.10.1`).
2. Install [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (highly recommended for full pipeline reproducibility).
3. Start running your analysis!

```bash
nextflow run nf_scrank/main.nf \
  --obj /path/to/your/seurat_object.rds \
  --column clone_annotation \
  --species human \
  --n_cells 3000 \
  --n_cores 32 \
  --target /path/to/targets.txt \
  --network genie3 \
  --outdir results \
  -profile singularity

``` 

### Inputs and References
NF_scRank requires two mandatory inputs:

Seurat Object (--obj): A fully processed .rds file containing normalized RNA assays and metadata annotations.

Target List (--target): A .txt file containing the list of genes of interest (one gene per line) for which the perturbation score will be calculated.

Note on cell numbers: If a specific cellular identity has fewer cells than the requested --n_cells parameter, the pipeline will automatically use all available cells for that identity.

### Outputs
If successfully run, the workflow will generate its primary output in the specified --outdir:

rank_scores/perbscore_all_targets.txt: A consolidated table containing the cell identity (cell_type), the evaluated gene (target), and its final regulatory importance (perb_score). This file is perfectly formatted to be loaded directly into R for visualization (e.g., ggplot2 barplots).

Other intermediate files (such as split matrices and raw GENIE3 weights) are temporarily stored in the work directory and can be retained or discarded based on standard Nextflow cache management.

### Credits
NF_scRank is developed and maintained by the Marques-Coelho Bioinformatics Lab(MCBLab).

### Citations
If you use this pipeline in your research, please cite:

https://www.sciencedirect.com/science/article/pii/S266637912400260X