FROM rocker/hadleyverse
MAINTAINER Stephen Eglen <sje30@cam.ac.uk>
RUN Rscript -e 'install.packages(c("splancs", "spatstat"))'
Run Rscript -e 'install.packages(c("sjedmin", "sjedrp", "sjevor","sjedist"), type="source", contriburl="http://damtp.cam.ac.uk/user/eglen/r/")'
