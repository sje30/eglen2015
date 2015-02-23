# eglen2015
Eglen 2015 review article.  I plan to put the R package here.


## Recompiling the paper

This R package depends on a few other packages, from CRAN and my
personal library.  The following sequence should install everything
you need:

    Rscript -e 'install.packages(c("splancs", "spatstat", "devtools", "knitr", "xtable"))'
    Rscript -e 'install.packages(c("sjedmin", "sjedrp", "sjevor","sjedist"), type="source", contriburl="http://damtp.cam.ac.uk/user/eglen/r/")'
    Rscript -e 'devtools::install_github("sje30/eglen2015",build_vignettes=TRUE)'
The last line should load this package.  Once it is installed, you can
then view the paper:

    vignette("eglen2015")

This does of course assume that your system already has R, latex, and
various unix tools.  That may not be the case; however, you can still
use the package through the docker system, see below.

### Docker

Once you have [docker](http://docker.com) installed on your system,
you can download and run this package using:

    docker run -d -p 8787:8787 sje30/eglen2015


Then visit the web page to start R (username and password are "rstudio"):

    http://localhost:8787/        ## linux
    http://192.168.59.103:8787/   ## mac, windows users


From that session, you can then view the vignette.

Thanks to the [Rocker](https://github.com/rocker-org) team for the
R-based docker images, on which this work is based.





