jmlr-paper.pdf: figure-all-cv.pdf
	rm -rf *.aux *.bbl
	pdflatex jmlr-paper
	bibtex jmlr-paper
	pdflatex jmlr-paper
	pdflatex jmlr-paper
figure-all-cv.pdf: figure-all-cv.R all.cv.RData
	R --vanilla < $<
all.cv.RData: all.cv.R all.modelSelection.RData
	R --vanilla < $<
jss-arxiv.pdf: jss-arxiv.Rnw
	rm -rf *.aux *.bbl
	R CMD Sweave jss-arxiv.Rnw
	pdflatex jss-arxiv
	bibtex jss-arxiv
	pdflatex jss-arxiv
	pdflatex jss-arxiv
	rm jss-arxiv.tex
jss-slides.pdf: jss-slides.tex jss-paper.pdf
	pdflatex jss-slides
jss-paper.pdf: jss-paper.Rnw jss-figure-more-likely-models-three-peaks.png jss-figure-target-intervals-models.pdf jss-figure-disk-memory-compare-speed.pdf jss-figure-data-peaks.tex jss-figure-label-error.pdf jss-figure-evaluations.tex jss-figure-variable-peaks.tex jss-refs.bib
	rm -rf *.aux *.bbl
	R CMD Sweave jss-paper.Rnw
	pdflatex jss-paper
	bibtex jss-paper
	pdflatex jss-paper
	pdflatex jss-paper
	rm jss-paper.tex
jss.evaluations.rds: jss.evaluations.R
	R --no-save < $<
jss.variable.peaks.rds: jss.variable.peaks.R
	R --no-save < $<
jss-figure-variable-peaks.tex: jss-figure-variable-peaks.R jss.variable.peaks.rds
	R --no-save < $<
jss-figure-evaluations.tex: jss-figure-evaluations.R jss.evaluations.rds
	R --no-save < $<
jss-figure-label-error.pdf: jss-figure-label-error.R
	R --no-save < $<
jss-figure-data-peaks.tex: jss-figure-data-peaks.R
	R --no-save < $<
jss.disk.memory.rds: jss.disk.memory.R
	R --vanilla < $<
jss-figure-disk-memory-compare-speed.pdf: jss-figure-disk-memory-compare-speed.R jss.disk.memory.rds
	R --no-save < $<
jss-figure-target-intervals-models.pdf: jss-figure-target-intervals-models.R target.intervals.models.csv
	R --no-save < $<
target.intervals.models.csv: target.intervals.models.R
	R --no-save < $<
jss-figure-more-likely-models-three-peaks.png: jss-figure-more-likely-models.R
	R --no-save < $<
HOCKING-RIGAILL-constrained-functional-pruning.pdf: HOCKING-RIGAILL-constrained-functional-pruning.tex refs.bib figure-1-min-operators.pdf figure-2-min-envelope.tex figure-PDPA-microbenchmark.pdf figure-PDPA-intervals.png figure-PDPA-timings.pdf figure-test-error-dots.pdf
	rm -rf *.aux *.bbl
	pdflatex HOCKING-RIGAILL-constrained-functional-pruning
	bibtex HOCKING-RIGAILL-constrained-functional-pruning
	pdflatex HOCKING-RIGAILL-constrained-functional-pruning
	pdflatex HOCKING-RIGAILL-constrained-functional-pruning
figure-large-margin/index.html: figure-large-margin.R
	R --no-save < $<
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
figure-Segmentor-PeakSeg.png: figure-Segmentor-PeakSeg.R
	R --no-save < $<
dp.peaks.NA.RData: dp.peaks.NA.R dp.peaks.matrices.RData
	R --no-save < $<
dp.peaks.matrices.RData: dp.peaks.matrices.R dp.peaks.error.RData PDPA.peaks.error.RData Segmentor.peaks.error.RData PDPA.infeasible.error.RData
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
figure-PDPA-timings.pdf: figure-PDPA-timings.R PDPA.timings.RData
	R --no-save < $<
PDPA.intervals.RData: PDPA.intervals.R
	R --no-save < $<
PDPA.model.check.RData: PDPA.model.check.R
	R --no-save < $<
PDPA.intervals.all.RData: PDPA.intervals.all.R
	R --no-save < $<
figure-PDPA-intervals.png: figure-PDPA-intervals.R PDPA.intervals.RData
	R --no-save < $<
figure-PDPA-intervals-all.png: figure-PDPA-intervals-all.R PDPA.intervals.all.RData
	R --no-save < $<
HOCKING-PeakSeg-functional-pruning-slides.pdf: HOCKING-PeakSeg-functional-pruning-slides.tex figure-macs-problem.png figure-min-train-error.pdf figure-min-undefined.pdf
	pdflatex $<
figure-macs-problem.png: figure-macs-problem.R
	R --no-save < $<
Segmentor.timings.RData: Segmentor.timings.R
	R --no-save < $<
Segmentor.peaks.error.RData: Segmentor.peaks.error.R
	R --no-save < $<
PDPA.peaks.error.RData: PDPA.peaks.error.R
	R --no-save < $<
PDPA.infeasible.error.RData: PDPA.infeasible.error.R PDPA.infeasible.RData
	R --no-save < $<
PDPA.infeasible.RData: PDPA.infeasible.R
	R --no-save < $<
figure-min-train-error.pdf: figure-min-train-error.R PDPA.peaks.error.RData Segmentor.peaks.error.RData
	R --no-save < $<
cosegData.timings.RData: cosegData.timings.R
	R --no-save < $<
figure-cosegData-timings.pdf: figure-cosegData-timings.R cosegData.timings.RData
	R --no-save < $<
figure-min-undefined.pdf: figure-min-undefined.R
	R --no-save < $<
macs.peaks.error.RData: macs.peaks.error.R
	R --no-save < $<
unsupervised.pdpa.RData: unsupervised.pdpa.R
	R --no-save < $<
unsupervised.RData: unsupervised.R
	R --no-save < $<
unsupervised.inf.RData: unsupervised.inf.R
	R --no-save < $<
test.error.RData: test.error.R unsupervised.RData unsupervised.inf.RData unsupervised.pdpa.RData dp.peaks.matrices.RData dp.peaks.sets.RData unsupervised.Segmentor.RData
	R --no-save < $<
dp.peaks.sets.RData: dp.peaks.sets.R dp.peaks.matrices.RData
	R --no-save < $<
figure-test-error-dots.pdf: figure-test-error-dots.R test.error.RData
	R --no-save < $<
unsupervised.Segmentor.RData: unsupervised.Segmentor.R
	R --no-save < $<
PDPA.cDPA.compare.RData: PDPA.cDPA.compare.R
	R --no-save < $<
figure-PDPA-cDPA-compare.pdf: figure-PDPA-cDPA-compare.R PDPA.cDPA.compare.RData
	R --no-save < $<
HOCKING-PeakSegFPOP-pipeline-slides.pdf: HOCKING-PeakSegFPOP-pipeline-slides.tex
	pdflatex HOCKING-PeakSegFPOP-pipeline-slides
PDPA.targets.RData: PDPA.targets.R PDPA.peaks.error.RData
	R --no-save < $<
problem.features.RData: problem.features.R
	R --no-save < $<
