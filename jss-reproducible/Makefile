FIGURES=jss-figure-more-likely-models-three-peaks.png jss-figure-target-intervals-models.pdf jss-figure-disk-memory-compare-speed.pdf jss-figure-data-peaks.tex jss-figure-label-error.pdf jss-figure-evaluations.tex jss-figure-variable-peaks.tex
jss-paper.pdf: jss-paper.Rnw jss-refs.bib ${FIGURES}
	rm -rf *.aux *.bbl
	R CMD Sweave jss-paper.Rnw
	pdflatex jss-paper
	bibtex jss-paper
	pdflatex jss-paper
	pdflatex jss-paper
	rm jss-paper.tex
jss.disk.memory.rds: jss.disk.memory.R
	R --vanilla < $<
jss.evaluations.rds: jss.evaluations.R
	R --vanilla < $<
jss.variable.peaks.rds: jss.variable.peaks.R
	R --vanilla < $<
jss.bench.models.csv: jss.bench.models.R 
	R --vanilla < $<
clean:
	rm ${FIGURES}
jss-figure-variable-peaks.tex: jss-figure-variable-peaks.R jss.variable.peaks.rds
	R --vanilla < $<
jss-figure-evaluations.tex: jss-figure-evaluations.R jss.evaluations.rds
	R --vanilla < $<
jss-figure-label-error.pdf: jss-figure-label-error.R
	R --vanilla < $<
jss-figure-data-peaks.tex: jss-figure-data-peaks.R
	R --vanilla < $<
jss-figure-disk-memory-compare-speed.pdf: jss-figure-disk-memory-compare-speed.R jss.disk.memory.rds
	R --vanilla < $<
jss-figure-target-intervals-models.pdf: jss-figure-target-intervals-models.R jss.bench.models.csv
	R --vanilla < $<
jss-figure-more-likely-models-three-peaks.png: jss-figure-more-likely-models.R
	R --vanilla < $<
