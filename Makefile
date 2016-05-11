HOCKING-PeakSegFPOP.pdf: HOCKING-PeakSegFPOP.tex refs.bib figure-unconstrained-PDPA-normal.pdf figure-unconstrained-FPOP-normal.pdf figure-constrained-PDPA-normal-grid.pdf figure-constrained-PDPA-normal-panels.pdf figure-less-more-min.tex figure-constrained-PDPA-normal-real.pdf figure-NA-timings.pdf
	pdflatex HOCKING-PeakSegFPOP
	bibtex HOCKING-PeakSegFPOP
	pdflatex HOCKING-PeakSegFPOP
	pdflatex HOCKING-PeakSegFPOP
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
