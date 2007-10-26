
FILES = RM.txt HOWTO.txt MOVE.txt
OBJS = $(FILES:.txt=.html)


all: $(OBJS) 

%.html: %.txt
	asciidoc $<

clean:
	rm -f *.html