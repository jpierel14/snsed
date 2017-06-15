# 2012.10.11 
# Steve Rodney

# This Makefile and the associated latex templates are useful
# for quick compilation of different versions of the same paper. 
# For example:
#
# make pdf : generate the compact version with emulateapj and pdflatex
# make manuscript : a manuscript version with lots of white space 
#     for co-authors to scribble comments
# make changetext : highlight in blue any text that you've flagged with
#      the \change{} command (e.g. for submitting a revised version back 
#      to the referee)
# make grayscale : use the grayscale version of figures, when available
# 
# make submission : the journal submission version, compiled with latex
#       (not pdflatex), using manuscript layout, eps figures, grayscale 
#       figs when available, red comments, and no blue text.


# Set the variable PAPER to the root name of your .tex document
PAPER= snsed

default: pdf

all: postscript manuscript grayscale changetext pdf
clean:
	rm -f *.aux *.log *.bbl *.blg *.lof *.lot *.toc *options.tex *~*

#Generate a compact pdf file using emulateapj
pdf: $(PAPER).tex
	pdflatex $(PAPER); \
	bibtex $(PAPER); \
	pdflatex $(PAPER); \
	pdflatex $(PAPER); \
	open $(PAPER).pdf


# quick update after minor changes
update: $(PAPER).tex
	pdflatex $(PAPER)

# re-open the file
open: $(PAPER).tex
	open $(PAPER).pdf

#Changes (marked with \change{} commands) are highlighted in blue text
changetext:
	echo "\\\changetexttrue" > $(PAPER)-options.tex; \
	echo "\\endinput" >> $(PAPER)-options.tex; \
	pdflatex $(PAPER); \
	bibtex $(PAPER); \
	pdflatex $(PAPER); \
	pdflatex $(PAPER); \
	rm $(PAPER)-options.tex; \
	mv $(PAPER).pdf $(PAPER)_changetext.pdf; \
	open $(PAPER)_changetext.pdf

