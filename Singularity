bootstrap:docker
From:rocker/hadleyverse

%post
# Stephen Eglen <sje30@cam.ac.uk>
apt-get update \

	&& apt-get install -t unstable -y --no-install-recommends
	texlive-bibtex-extra
Rscript -e 'install.packages(c("splancs", "spatstat", "knitr",
"xtable"))'
Rscript -e 'install.packages(c("sjedmin", "sjedrp",
"sjevor","sjedist"), type="source",
contriburl="http://damtp.cam.ac.uk/user/eglen/r/")'
Rscript -e
'devtools::install_github("sje30/eglen2015",build_vignettes=TRUE)'
