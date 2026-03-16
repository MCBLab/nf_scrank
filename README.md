# **NF_scRank** 🧬🥇
This repository is intended to document the development of a benchmark based on pipeline integration (Nextflow) and a Single-Cell expression ranking algorithm (scRank).

## **Script usage**

### Prerequisites 
Before running the pipeline, ensure you have the following installed in your environment:
* [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) 
* [Singularity](https://sylabs.io/guides/3.0/user-guide/index.html)

### Usage
The pipeline requires a processed Seurat object (`.rds`) and a target genes list (`.txt`). It automatically handles downsampling, split processes and scoring consolidation.

Here is a basic execution example using Singularity:

```bash
nextflow run MCBLab/nf_scrank \
  --obj /path/to/your/seurat_object.rds \
  --column clone_annotation \
  --species human \
  --n_cells 3000 \
  --n_cores 32 \
  --target /path/to/targets.txt \
  --outdir results \
  -profile singularity \
  -resume

### Key Parameters

### Output
The pipeline generates a consolidated .txt file with the perturbation scores for each target gene across all specified cell identities, ready for downstream visualization in R


###references

https://www.sciencedirect.com/science/article/pii/S266637912400260X