LATEX=pdflatex
BIB=biber
RM=rm -f -r

all: report | clean_not_pdf # build, then clean up

report: report.tex
	$(LATEX) report.tex
	$(BIB) report
	$(LATEX) report.tex
	$(LATEX) report.tex	

clean_not_pdf:
	$(RM) prv_report.* _region_.* *.dvi *.log *.aux *.run.xml *.bcf *.blg *.bbl *.out

clean: clean_not_pdf
	$(RM) *.pdf
