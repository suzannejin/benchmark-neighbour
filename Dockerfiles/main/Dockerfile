FROM r-base:4.1.3

# restore r environment
RUN git clone https://github.com/suzannejin/benchmark-neighbour.git
RUN R -e "renv::restore()"
RUN R -e "reticulate::virtualenv_install('normalisr_python_env', packages = 'normalisr')"
