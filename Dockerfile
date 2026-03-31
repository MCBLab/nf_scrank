FROM satijalab/seurat:5.4.0

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libglpk-dev \
    libbz2-dev \
    liblzma-dev \
    libzstd-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('igraph', 'remotes', 'qs', 'harmony', 'fastcluster'), repos='http://cran.us.r-project.org')"

RUN R -e "BiocManager::install(c( \
    'WGCNA', \
    'impute', \
    'preprocessCore', \
    'GO.db', \
    'AnnotationDbi', \
    'qvalue', \
    'ComplexHeatmap', \
    'GeneOverlap', \
    'UCell' \
), ask=FALSE, update=FALSE)"

RUN R -e "options(timeout = 600); remotes::install_github('smorabit/hdWGCNA', ref='dev', upgrade='never'); if (!requireNamespace('hdWGCNA', quietly = TRUE)) stop('ERRO CRITICO: hdWGCNA nao instalou!')"

ENV R_MAX_VSIZE=16Gb
WORKDIR /workspace

CMD ["R"]