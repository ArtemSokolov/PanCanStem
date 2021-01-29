FROM rocker/tidyverse:4.0.1

RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install("biomaRt")'
RUN R -e 'install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))'
RUN R -e 'devtools::install_github("ArtemSokolov/gelnet")'

COPY ./*.R /app/
