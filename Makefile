HOCKING-PeakSegFPOP.pdf: HOCKING-PeakSegFPOP.tex refs.bib figure-unconstrained-PDPA-normal.pdf figure-unconstrained-FPOP-normal.pdf
	pdflatex HOCKING-PeakSegFPOP
	bibtex HOCKING-PeakSegFPOP
	pdflatex HOCKING-PeakSegFPOP
	pdflatex HOCKING-PeakSegFPOP
figure-unconstrained-PDPA-normal.pdf: figure-unconstrained-PDPA-normal.R
	R --no-save < $<
figure-unconstrained-FPOP-normal.pdf: figure-unconstrained-FPOP-normal.R
	R --no-save < $<

