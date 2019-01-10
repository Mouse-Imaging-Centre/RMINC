all:
	Rscript -e "rmarkdown::render('index.Rmd')"

.PHONY: all
