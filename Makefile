install:
	Rscript -e 'devtools::install(build_vignettes=TRUE)'

dockerfile-lint:
	docker run --rm -i hadolint/hadolint < Dockerfile
