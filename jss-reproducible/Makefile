FIGURES=jss-figure-more-likely-models-three-peaks.png jss-figure-target-intervals-models.pdf jss-figure-disk-memory-compare-speed.pdf jss-figure-data-peaks.pdf jss-figure-label-error.pdf jss-figure-evaluations.pdf jss-figure-variable-peaks.pdf
jss-figures/DONE: ${FIGURES}
	mkdir -p jss-figures
	cp jss-figure-evaluations-concave.pdf jss-figures/figure-3-left.pdf
	cp jss-figure-evaluations-concave-zoom.pdf jss-figures/figure-3-right.pdf
	cp jss-figure-more-likely-models-three-peaks.png jss-figures/figure-4-left.png
	cp jss-figure-more-likely-models-three-peaks-zoom.png jss-figures/figure-4-right.png
	cp jss-figure-target-intervals-models.pdf jss-figures/figure-5-left.pdf
	cp jss-figure-target-intervals-models-computation.pdf jss-figures/figure-5-right.pdf
	cp jss-figure-disk-memory-compare-speed.pdf jss-figures/figure-6-left.pdf
	cp jss-figure-disk-memory-compare-speed-penalty.pdf jss-figures/figure-6-right.pdf
	cp jss-figure-variable-peaks.pdf jss-figures/figure-7-left.pdf
	cp jss-figure-variable-peaks-zoom.pdf jss-figures/figure-7-right.pdf
	cp jss-figure-label-error.pdf jss-figures/figure-8.pdf
	cp jss-figure-data-peaks.pdf jss-figures/figure-9.pdf
	cp jss-figure-evaluations.pdf jss-figures/figure-10-left.pdf
	cp jss-figure-evaluations-computation.pdf jss-figures/figure-10-right.pdf
	touch jss-figures/DONE 
jss.disk.memory.rds: jss.disk.memory.R
	R --vanilla < $<
jss.evaluations.rds: jss.evaluations.R
	R --vanilla < $<
jss.variable.peaks.rds: jss.variable.peaks.R
	R --vanilla < $<
jss.bench.models.csv: jss.bench.models.R 
	R --vanilla < $<
clean:
	rm -rf *.tex *.pdf *.png library *~ logs Mono27ac tikzMetricsJSS monitor.sh auto jss-figures
jss-figure-variable-peaks.pdf: jss-figure-variable-peaks.R jss.variable.peaks.rds
	R --vanilla < $<
jss-figure-evaluations.pdf: jss-figure-evaluations.R jss.evaluations.rds
	R --vanilla < $<
jss-figure-label-error.pdf: jss-figure-label-error.R
	R --vanilla < $<
jss-figure-data-peaks.pdf: jss-figure-data-peaks.R
	R --vanilla < $<
jss-figure-disk-memory-compare-speed.pdf: jss-figure-disk-memory-compare-speed.R jss.disk.memory.rds
	R --vanilla < $<
jss-figure-target-intervals-models.pdf: jss-figure-target-intervals-models.R jss.bench.models.csv
	R --vanilla < $<
jss-figure-more-likely-models-three-peaks.png: jss-figure-more-likely-models.R
	R --vanilla < $<
