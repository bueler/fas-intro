all: fas.pdf

publish.zip: fas.pdf fas.bbl
	rm -rf publish/
	mkdir publish/
	mkdir publish/figs/
	cp -f figs/*.pdf publish/figs/
	mkdir publish/tikz/
	cp -f tikz/*.tex publish/tikz/
	cp -f fas.bib fas.bbl fas.tex publish/
	zip -r publish publish

_figures := onehat hatcombination converge optimal tme
figures := $(addprefix figs/, $(addsuffix .pdf, $(_figures)))

fas.pdf: fas.tex fas.bib $(figures)
	pdflatex fas
	bibtex fas
	pdflatex fas
	pdflatex fas

$(figures): figs/hatfig.py figs/converge.py figs/optimal.py figs/tme.py
	(cd figs/ && python3 hatfig.py && python3 converge.py && python3 optimal.py && python3 tme.py)
	for X in $(figures); do \
	     (pdfcrop $$X $$X); \
	done

.PHONY: clean

clean:
	@rm -f *~ *.out *.aux *.log *.blg *.snm *.toc *.nav *.vrb *.bbl *.synctex.gz
	@rm -f $(figures)

