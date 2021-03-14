# adipo_deconv
Code for the preprint "A small fraction of progenitors differentiate into mature adipocytes due to constraints on the cell structure change"

## Setting up the docker environment

The analysis was run on a [docker](https://hub.docker.com/r/bcmslab/adipo_deconv/)
image based on the the latest **rocker/verse**.
Other R packages were added to the image and were made available as an image 
that can be obtained and launched on any local machine running
[docker](https://hub.docker.com/r/bcmslab/adipo_deconv/).

```bash
$ docker pull bcmslab/adipo_deconv:latest
$ docker run -it bcmslab/adipo_deconv:latest bash
```

## Obtaining the source code

The source code is hosted publicly on this repository in the form of a research
compendium. This includes the scripts to reproduce the figures and tables in 
this manuscript.

```bash
$ git clone https://github.com/BCMSLab/adipo_deconv
```

## Runing the analysis

In the directory `adipo_deconv`, 
- Run `scripts/initial_analysis.R` to reproduce the analysis
- Knit the `adipo_deconv.Rmd` to reproduce the figures and tables in the manuscript

```bash
$ cd adipo_deconv
$ Rscript scripts/initial_analysis.R
$ Rscript -e "knit::knit('adipo_deconv.Rmd')"
```

## Details of the R environment
The version of **R** that was used to perform this analysis is the 4.02
(2020-06-22) on `x86\_64-pc-linux-gnu`.

## More

This manuscript was released as a preprint under the title [A small fraction of progenitors differentiate into mature adipocytes due to constraints on the cell structure change](https://doi.org/10.1101/2021.03.10.434887)
