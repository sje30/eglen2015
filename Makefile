install:
	Rscript -e 'devtools::install(build_vignettes=TRUE)'
	##cd vignettes; Rscript -e "knitr::knit2pdf('eglen2015.Rnw')"

