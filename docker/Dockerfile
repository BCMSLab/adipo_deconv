FROM bioconductor/bioconductor_docker:latest

RUN sudo apt-get -y update
RUN sudo apt-get -y install libproj-dev

RUN Rscript -e "install.packages('tidyverse')"
RUN Rscript -e "install.packages('reshape2')"
RUN Rscript -e "install.packages('pdist')"
RUN Rscript -e "install.packages('circlize')"
RUN Rscript -e "install.packages('xtable')"
RUN Rscript -e "install.packages('deSolve')"
RUN Rscript -e "install.packages('scales')"
RUN Rscript -e "install.packages('cowplot')"
RUN Rscript -e "install.packages('ggalt')"
RUN Rscript -e "install.packages('gridExtra')"
RUN Rscript -e "install.packages('png')"
RUN Rscript -e "install.packages('BiocManager')"

RUN Rscript -e "options(repos = BiocManager::repositories())"

RUN Rscript -e "BiocManager::install('SummarizedExperiment')"
RUN Rscript -e "BiocManager::install('org.Mm.eg.db')"
RUN Rscript -e "BiocManager::install('GO.db')"
RUN Rscript -e "BiocManager::install('ComplexHeatmap')"
RUN Rscript -e "BiocManager::install('clusterProfiler')"
RUN Rscript -e "BiocManager::install('limma')"
RUN Rscript -e "BiocManager::install('PSEA')"
RUN Rscript -e "BiocManager::install('DESeq2')"
RUN Rscript -e "BiocManager::install('debCAM')"
RUN Rscript -e "install.packages('usedist')"
RUN Rscript -e "install.packages('Hmisc')"
RUN Rscript -e "install.packages('pcr')"
RUN Rscript -e "install.packages('GGally')"
