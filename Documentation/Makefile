all: Manual.ps

Manual.tex : part1.tex part2.tex part3.tex references.bib 
	touch Manual.tex

Manual.dvi : Manual.tex
	latex Manual.tex
	bibtex Manual
	makeindex aedl.idx
	makeindex norm.idx
	makeindex api.idx
	makeindex constants.idx
	remove_spaces constants.ind
	remove_spaces api.ind
	remove_spaces norm.ind
	remove_spaces aedl.ind
	latex Manual.tex
	latex Manual.tex

Manual.ps : Manual.dvi
	dvips Manual.dvi

.PHONY: clean

clean:
	rm -f *.aux *.log *.ind *.bbl *.ilg *.blg Manual.dvi Manual.ps *.toc 

