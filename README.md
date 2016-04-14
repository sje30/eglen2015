# Introduction

[![Build Status](https://travis-ci.org/sje30/eglen2015.svg?branch=master)](https://travis-ci.org/sje30/eglen2015)



This is the home page for the following review:

Eglen SJ (2016) Bivariate spatial point patterns in the retina: a
reproducible review. Journal de la Société Française de Statistique
157:33–48.
[PDF](http://journal-sfds.fr/index.php/J-SFdS/article/view/518/490)

This package contains all the material needed to regenerate the
article for itself.  Some key parts of the package are:

* [vignettes/eglen2015.Rnw](vignettes/eglen2015.Rnw): the source file
  for the article in LaTeX format.
* [inst/extdata](inst/extdata): a folder containing all the data files
studied in this article.



## Recompiling the paper

This R package depends on a few other packages, from CRAN and my
personal library.  The following sequence should install everything
you need:

    Rscript -e 'install.packages(c("splancs", "spatstat", "devtools", "knitr", "xtable"))'
    Rscript -e 'install.packages(c("sjedmin", "sjedrp", "sjevor","sjedist"), type="source", contriburl="http://damtp.cam.ac.uk/user/eglen/r/")'
    Rscript -e 'devtools::install_github("sje30/eglen2015",build_vignettes=TRUE)'
The last line should load this package.  Once it is installed, you can
then view the paper, or view the knitr document that created the paper:

    vignette("eglen2015")
	eglen2015:::edit()
	
This does of course assume that your system already has R, latex, and
various unix tools.  That may not be the case; however, you can still
use the package through the Docker system, see next.



### Docker

Once you have [docker](http://docker.com) installed on your system,
you can download and run this package using:

    docker run -d -p 8787:8787 sje30/eglen2015

(View [sje30/eglen2015](https://registry.hub.docker.com/u/sje30/eglen2015/)
to check the status of this Docker package.)

Then visit the web page to start R (username and password are "rstudio"):

    http://localhost:8787/        ## linux
    http://192.168.99.100:8787/   ## mac, windows users

The IP address for mac/windows may vary; you can check it by running
the command:

	docker-machine ip default

In the File pane (bottom right), you can click on the file `README.R`
to get started.  In the R Console (bottom left), you can view the
vignette (in PDF form) with `vignette("eglen2015")`.  You can edit the
source of the paper with the command `eglen2015:::edit()` or re-build
it from source with the command `tools::buildVignettes("eglen2015")`.
Thanks to the [Rocker](https://github.com/rocker-org) team for the
R-based docker images, on which this work is based.



