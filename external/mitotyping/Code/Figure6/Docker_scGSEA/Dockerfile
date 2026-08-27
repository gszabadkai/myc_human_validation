# Use RStudio-enabled image with R 4.3.0
FROM rocker/rstudio:4.3.0

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install Linux system dependencies (required for many R packages)
RUN apt-get update && apt-get install -y \
    libgsl-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    #libsqlite3-dev \
    #libpq-dev \
    #libfontconfig1-dev \
    #libfreetype6-dev \
    build-essential \
    git \
    libblas-dev \
    liblapack-dev \
    gfortran \
    libffi-dev \
    libpng-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libpcre2-dev \
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    libpython3-dev \
    wget \
    curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


# Upgrade pip
RUN python3 -m pip install --upgrade pip setuptools wheel

# Install Python packages for reticulate and leiden
RUN pip3 install igraph leidenalg pandas umap-learn

# Set R default options (libcurl + CRAN mirror)
RUN echo 'options(download.file.method = "libcurl")' >> /usr/local/lib/R/etc/Rprofile.site \
 && echo 'options(repos = c(CRAN = "https://cloud.r-project.org"))' >> /usr/local/lib/R/etc/Rprofile.site

# Set C++14 as default for R packages
#RUN echo "CXX = g++ -std=c++14" >> /usr/local/lib/R/etc/Makeconf \
# && echo "CXX14 = g++ -std=c++14" >> /usr/local/lib/R/etc/Makeconf \
# && echo "CXX14STD = -std=c++14" >> /usr/local/lib/R/etc/Makeconf \
# && echo "CXX11 = g++ -std=c++14" >> /usr/local/lib/R/etc/Makeconf

# Fix permissions so rstudio user can install packages
RUN chown -R rstudio:rstudio /usr/local/lib/R /usr/local/lib/R/site-library

# Install R package management helpers
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org')" \
 && R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org')"

# Install core Bioconductor packages
RUN R -e "BiocManager::install(c('AnnotationDbi', 'AnnotationHub'), dependencies = TRUE)"

# Install CRAN packages
RUN R -e "install.packages(c( \
    'uwot', 'igraph', 'ggplot2', 'irlba', 'ggrepel', 'Rtsne', \
    'fastmatch', 'RSpectra', 'Matrix', 'babelgene', 'reshape2', \
    'RcppParallel', 'Rcpp', 'MASS', 'locfit', 'leiden', \
    'KernelKnn', 'mgcv', 'pointr', 'knitr', 'utils', 'BH'), repos='https://cloud.r-project.org')"

# Install remaining Bioconductor packages
RUN R -e "BiocManager::install(c('edgeR', 'GSVA', 'limma', 'harmony', 'BiocParallel', 'ensembldb', 'sva', 'msigdbr'), ask = FALSE)"

# Install reticulate
RUN R -e "install.packages('reticulate', repos='https://cloud.r-project.org')"

# Install GitHub packages
RUN R -e "remotes::install_github('zdebruine/RcppML')" \
 && R -e "remotes::install_github('alserglab/fgsea')" \
 && R -e "remotes::install_github('gambalab/gficf')"

# Set working directory
WORKDIR /home/rstudio

# Default user
RUN chown -R rstudio:rstudio /home/rstudio