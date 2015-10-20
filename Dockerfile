FROM rocker/hadleyverse
MAINTAINER Stephen Eglen <sje30@cam.ac.uk>
RUN apt-get update \
	&& apt-get install -t unstable -y --no-install-recommends texlive-bibtex-extra
RUN Rscript -e 'install.packages(c("splancs", "spatstat", "knitr", "xtable"))'
RUN Rscript -e 'install.packages(c("sjedmin", "sjedrp", "sjevor","sjedist"), type="source", contriburl="http://damtp.cam.ac.uk/user/eglen/r/")'
RUN Rscript -e 'devtools::install_github("sje30/eglen2015",build_vignettes=TRUE)'
RUN Rscript -e 'file.copy(system.file("doc/README.R", package="eglen2015"), "/home/rstudio/README.R")'

## To rebuild:
## docker build -t sje30/eglen2015 https://raw.githubusercontent.com/sje30/eglen2015/master/Dockerfile

## texlive-bibtex-extra is required for the breakcites.sty package
## which in turn is needed by the JSFDS package.
