install:
	Rscript -e 'devtools::install(build_vignettes=FALSE)'
	cd vignettes; Rscript -e "knitr::knit2pdf('eglen2015.Rnw')"

