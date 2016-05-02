HOCKING-PeakSegFPOP.pdf: HOCKING-PeakSegFPOP.tex refs.bib figure-unconstrained-PDPA-normal.pdf figure-unconstrained-FPOP-normal.pdf figure-constrained-PDPA-normal-grid.pdf figure-constrained-PDPA-normal.pdf
	pdflatex HOCKING-PeakSegFPOP
	bibtex HOCKING-PeakSegFPOP
	pdflatex HOCKING-PeakSegFPOP
	pdflatex HOCKING-PeakSegFPOP
figure-constrained-PDPA-normal-grid.pdf: figure-constrained-PDPA-normal-grid.R
	R --no-save < $<
figure-constrained-PDPA-normal.pdf: figure-constrained-PDPA-normal.R
	R --no-save < $<
figure-unconstrained-PDPA-normal.pdf: figure-unconstrained-PDPA-normal.R
	R --no-save < $<
figure-unconstrained-FPOP-normal.pdf: figure-unconstrained-FPOP-normal.R
	R --no-save < $<

