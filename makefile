
all: landing walkthroughs

clean:
	if [ -d docs ]; then rm -r docs; fi

landing:
	mkdir -p docs
	R -e "rmarkdown::render('index.md', output_dir = 'docs')"
	R -e "rmarkdown::render('analysis_template.Rmd', output_dir = 'docs')"
	R -e "rmarkdown::render('code_of_conduct.md', output_dir = 'docs')"

walkthroughs: 
	mkdir -p docs/tutorials
	R -e "f <- list.files('tutorials', '.Rmd', full.name = TRUE); lapply(f, rmarkdown::render, output_dir = 'docs/tutorials')"


