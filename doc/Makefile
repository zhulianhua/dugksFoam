VERS = 1.1
FILES = biblio.tex installation.tex demo.tex contrib.tex lshort.sty mylayout.sty title.tex \
	dugksFoam.tex overview.tex fancyhea.sty changelog.tex\
	lshort-base.tex appendix.tex

# Define some variables
PDFLATEX=pdflatex
MAKEINDEX=makeindex

# The default targets
all: dugksFoam.pdf 

dugksFoam.pdf: $(FILES)
	-mkdir pdfbuild
	(T1FONTS=.:`pwd`/eurofont: && export T1FONTS && TEXINPUTS=.:`pwd`:`pwd`/euro:`pwd`/oberdiek:${TEXINPUTS:-:}&&export TEXINPUTS&& cd pdfbuild&& \
	$(PDFLATEX) dugksFoam&& $(PDFLATEX) dugksFoam&& \
	#$(MAKEINDEX) -s dugksFoam.ist dugksFoam&&$(PDFLATEX) dugksFoam&& \
	(thumbpdf --resolution 10 dugksFoam.pdf && $(PDFLATEX) dugksFoam)&& \
	mv dugksFoam.pdf .. )
	rm pdfbuild/*

clean:
	rm -rf texbuild pdfbuild
