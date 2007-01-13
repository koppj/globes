# This Makefile produces the Manual for GLoBES
#
# 2004 Patrick Huber
#
# It takes care of the multiple indices and calls the script
# remove_spaces to fix a kind of font problem with the makeindex
# command.
#
# It also produces a pdf file on request (make pdf) by using
# pdflatex. This requires to have all eps files converted to 
# pdf file, which is done by epstopdf. If you add any new eps
# files, they automatically will be converted to pdf files.
# In order to use the same tex sources for both ps and pdf output
# please use
# \includegraphics{filename}
# *without* the .eps or .pdf extension.
# The resulting pdf file should be of high quality and well readable on
# screen.

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
	./pdfmanual

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
