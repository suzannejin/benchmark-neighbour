FROM r-base:4.1.3

RUN apt-get update \
    && apt-get install -y git libcurl4-openssl-dev libssl-dev libxml2-dev libgeos-dev cmake \
    && git clone https://github.com/suzannejin/benchmark-neighbour.git \
    && cd benchmark-neighbour/benchmark \
    && R -e "renv::restore()" \
    && R -e "reticulate::install_miniconda()" \
    && R -e "reticulate::virtualenv_install('normalisr_python_env', packages = 'normalisr')"
