FROM satijalab/seurat:5.4.0

# Install Bioconductor dependencies
RUN R -e "BiocManager::install(c( \
    'BiocGenerics', \
    'SummarizedExperiment' \
), ask=FALSE)"

# Install scTenifoldNet from GitHub
RUN R -e "devtools::install_github('cailab-tamu/scTenifoldNet')"

# Optional: set parallelization defaults
ENV R_MAX_VSIZE=8Gb

WORKDIR /workspace

CMD ["R"]