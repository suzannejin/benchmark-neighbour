FROM r-base:4.1.3

WORKDIR /opt

# Install libraries required for R packages
RUN apt update && \
    apt install -y libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfftw3-dev \
    libtiff-dev \
    libgsl-dev \
    libgeos-dev \
    cmake \
    git

# Install renv
ENV RENV_VERSION 0.15.4
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" && \
    R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

# Add ENV variables for Renv cache (not symlink)
RUN mkdir /opt/renvcache && \
 echo RENV_PATHS_CACHE='/opt/renvcache' >> $(R RHOME)/etc/Renviron.site

# set the project directory to a variable
ENV RENV_DIR '/opt/benchmark-neighbour/benchmark'

# restore R packages and move the symlinks
RUN git clone https://github.com/suzannejin/benchmark-neighbour.git && \
    cd ${RENV_DIR} && \
    R -e "renv::restore()" && \
    R -e "reticulate::install_miniconda()" && \
    R -e "reticulate::virtualenv_install('normalisr_python_env', packages = 'normalisr')" && \
    mv ${RENV_DIR}/renv/library/R-4.1/x86_64-pc-linux-gnu /opt/Rlibsymlinks && \
    echo "R_LIBS=/opt/Rlibsymlinks" >> $(R RHOME)/etc/Renviron.site && \
    rm -rf /opt/benchmark-neighbour
