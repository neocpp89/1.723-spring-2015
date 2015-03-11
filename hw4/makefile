TEXSRC = writeup.tex
PDF = $(TEXSRC:.tex=.pdf)

.PHONY: clean

all: $(PDF)

clean:
	latexmk -C $(TEXSRC);

$(PDF): $(TEXSRC) $(SECSRC) $(BIBSRC)
	latexmk -pdf -pdflatex="pdflatex --shell-escape %O %S" $(TEXSRC); latexmk -c $(TEXSRC);

