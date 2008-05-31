# Makefile fuer A0-Poster

POSTERNAME=globes-poster

.PHONY: clean

psA4: ps
	pstops -h3402 -w2594 "1:0@0.25(0cm,0cm)" ${POSTERNAME}.ps ${POSTERNAME}-a4.ps

ps: dvi
	dvips ${POSTERNAME}.dvi

pdf:
	pdflatex ${POSTERNAME}.tex

dvi: ${POSTERNAME}.tex
	latex ${POSTERNAME}.tex

clean:
	rm -f *.ps *.pdf *.aux *.log *.dvi


