all: Manual.ps

pdfpics := $(patsubst %.eps,%.pdf,$(wildcard *.eps))

$(pdfpics) : %.pdf : %.eps
	epstopdf $< --outfile=$@ 

epspdf : $(pdfpics)

Manual.tex : part1.tex part2.tex part3.tex references.bib appendix.tex
	touch Manual.tex

Manual.dvi : Manual.tex
	latex Manual.tex
	bibtex Manual
	makeindex aedl.idx
	makeindex norm.idx
	makeindex api.idx
	makeindex constants.idx
	./remove_spaces constants.ind
	./remove_spaces api.ind
	./remove_spaces norm.ind
	./remove_spaces aedl.ind
	latex Manual.tex
	latex Manual.tex

Manual.pdf : Manual.dvi epspdf
	pdfmanual

pdf: Manual.pdf

Manual.ps : Manual.dvi
	dvips Manual -o Manual.ps

.PHONY: clean fast epspdf pdf

clean:
	rm -f *.aux *.log *.ind *.bbl *.ilg *.blg Manual.dvi Manual.ps *.toc\
Manual.pdf *.idx

fast:
	latex Manual
	dvips Manual -o Manual.ps
