HOCKING-RIGAILL-constrained-functional-pruning.pdf: HOCKING-RIGAILL-constrained-functional-pruning.tex refs.bib figure-1-min-operators.pdf figure-2-min-envelope.tex figure-PDPA-microbenchmark.pdf figure-PDPA-intervals.pdf
	rm -rf *.aux *.bbl
	pdflatex HOCKING-RIGAILL-constrained-functional-pruning
	bibtex HOCKING-RIGAILL-constrained-functional-pruning
	pdflatex HOCKING-RIGAILL-constrained-functional-pruning
	pdflatex HOCKING-RIGAILL-constrained-functional-pruning
figure-2-min-envelope.tex: figure-2-min-envelope.R
	R --no-save < $<
figure-1-min-operators.pdf: figure-1-min-operators.R
	R --no-save < $<
figure-cDPA-PDPA-all/index.html: figure-cDPA-PDPA-all.R
	R --no-save < $<
figure-cDPA-PDPA.pdf: figure-cDPA-PDPA.R
	R --no-save < $<
HOCKING-notes.pdf: HOCKING-notes.tex refs.bib figure-unconstrained-PDPA-normal.pdf figure-unconstrained-FPOP-normal.pdf figure-constrained-PDPA-normal-grid.pdf figure-constrained-PDPA-normal-panels.pdf figure-less-more-min.tex figure-constrained-PDPA-normal-real.pdf figure-NA-timings.pdf
	pdflatex HOCKING-notes
	bibtex HOCKING-notes
	pdflatex HOCKING-notes
	pdflatex HOCKING-notes
figure-constrained-PDPA-normal-real.pdf: figure-constrained-PDPA-normal-real.R
	R --no-save < $<
figure-less-more-min.tex: figure-less-more-min.R
	R --no-save < $<
figure-constrained-PDPA-normal-panels.pdf: figure-constrained-PDPA-normal-panels.R
	R --no-save < $<
figure-constrained-PDPA-normal-grid.pdf: figure-constrained-PDPA-normal-grid.R
	R --no-save < $<
figure-constrained-PDPA-normal.pdf: figure-constrained-PDPA-normal.R
	R --no-save < $<
figure-unconstrained-PDPA-normal.pdf: figure-unconstrained-PDPA-normal.R
	R --no-save < $<
figure-unconstrained-FPOP-normal.pdf: figure-unconstrained-FPOP-normal.R
	R --no-save < $<
figure-NA-timings.pdf: figure-NA-timings.R dp.peaks.NA.RData
	R --no-save < $<
## Copied from PeakSeg paper.
dp.peaks.NA.RData: dp.peaks.NA.R dp.peaks.matrices.RData
	R --no-save < $<
dp.peaks.matrices.RData: dp.peaks.matrices.R dp.peaks.error.RData
	R --no-save < $<
dp.peaks.error.RData: dp.peaks.error.R dp.peaks.RData
	R --no-save < $<
dp.peaks.RData: dp.peaks.R dp.timings.RData
	R --no-save < $<
dp.timings.RData: dp.timings.R
	R --no-save < $<
dp.timings.reverse.RData: dp.timings.reverse.R
	R --no-save < $<
PDPA.timings.RData: PDPA.timings.R
	R --no-save < $<
PDPA.microbenchmark.RData: PDPA.microbenchmark.R
	R --no-save < $<
figure-PDPA-microbenchmark.pdf: figure-PDPA-microbenchmark.R PDPA.microbenchmark.RData
	R --no-save < $<
PDPA.intervals.RData: PDPA.intervals.R
	R --no-save < $<
figure-PDPA-intervals.pdf: figure-PDPA-intervals.R PDPA.intervals.RData
	R --no-save < $<
