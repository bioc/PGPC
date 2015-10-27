all:
	$(R_HOME)/bin/R CMD pdflatex PGPC.tex
	$(R_HOME)/bin/R CMD bibtex PGPC
	$(R_HOME)/bin/R CMD pdflatex PGPC.tex
	$(R_HOME)/bin/R CMD pdflatex PGPC.tex
	make clean

clean:
	make texclean
	rm -rf result
	rm -rf cache
	rm -rf figure
	rm -rf src/result

texclean:
	rm -f *.tex
	rm -f *.aux
	rm -f *.log
	rm -f *.out
	rm -f *.toc
	rm -f *.bbl
	rm -f *.blg
