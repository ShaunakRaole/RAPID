FROM rocker/r-ver:4.3.3

LABEL maintainer="Shaunak Raole <sraole3@gatech.edu>"
LABEL description="RAPID — Reproducible Analysis Pipeline for Intensity-based proteomics Data"
LABEL version="1.0"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    pandoc \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# CRAN packages — pinned snapshot for reproducibility
RUN Rscript -e "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/2024-06-01')); \
    install.packages(c( \
        'data.table', \
        'ggplot2', \
        'ggrepel', \
        'rmarkdown', \
        'knitr', \
        'kableExtra', \
        'bit64' \
    ), dependencies = TRUE)"

# Bioconductor packages
RUN Rscript -e "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/2024-06-01')); \
    if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
    BiocManager::install(version = '3.18', ask = FALSE, update = FALSE); \
    BiocManager::install(c( \
        'limma', \
        'DEqMS', \
        'clusterProfiler', \
        'org.Hs.eg.db', \
        'enrichplot' \
    ), ask = FALSE, update = FALSE)"

# Verify all packages installed
RUN Rscript -e " \
    pkgs <- c('data.table', 'ggplot2', 'ggrepel', 'rmarkdown', \
              'knitr', 'kableExtra', 'bit64', \
              'limma', 'DEqMS', 'clusterProfiler', \
              'org.Hs.eg.db', 'enrichplot'); \
    for (p in pkgs) { \
        if (!requireNamespace(p, quietly = TRUE)) \
            stop(paste('Package not installed:', p)); \
        cat(sprintf('OK: %s %s\n', p, as.character(packageVersion(p)))) \
    }"

WORKDIR /pipeline
CMD ["Rscript", "--version"]
