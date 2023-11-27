.PHONY: programs 

BINROOT = ..
RIGID2D=$(BINROOT)/Examples/MultiResImageRegistration2D
GRADMAP=$(BINROOT)/Examples/NGFGradientMap
Image2Grad=$(BINROOT)/Examples/ImageToNGF

NAME = NGFImageMetric

BIBTEXDATABSENAME = InsightJournal

GRADIMAGES= fscalar.png fcross.png fscdelta.png \
	rscalar.png rcross.png rscdelta.png \
	self-scdelta.png self-scalar.png self-cross.png \
	algscdelta1-2.png algovscdelta1-2.png ovscdelta1-2.png \
	algscalar1-2.png algovscalar1-2.png ovscalar1-2.png\
	algcross1-2.png algovcross1-2.png ovcross1-2.png \
	algscdelta2-1.png algovscdelta2-1.png ovscdelta2-1.png \
	algscalar2-1.png algovscalar2-1.png ovscalar2-1.png\
	algcross2-1.png algovcross2-1.png ovcross2-1.png \
	src1grad.png src2grad.png 

all: pdf 



pdf: ${NAME}.pdf


bbl: ${NAME}.bbl


${NAME}.pdf: ${NAME}.tex ${NAME}.bbl $(GRADIMAGES)
	pdflatex ${NAME}.tex
	pdflatex ${NAME}.tex


${NAME}.bbl: ${BIBTEXDATABSENAME}.bib ${NAME}.aux
	bibtex ${NAME}


${NAME}.aux: ${NAME}.tex  ${BIBTEXDATABSENAME}.bib $(GRADIMAGES)
	pdflatex ${NAME}.tex


clean:
	rm ${NAME}.pdf ${NAME}.aux ${NAME}.log ${NAME}.bbl ${NAME}.blg ${NAME}.toc ${NAME}.out  ${NAME}.brf

r%.png: src1.png src2.png $(GRADMAP)
	$(GRADMAP) src2.png src1.png $(shell echo $@|sed -s "s/r\([a-z]*\).png/\1/") $@ 


f%.png: src1.png src2.png $(GRADMAP)
	$(GRADMAP) src1.png src2.png $(shell echo $@|sed -s "s/f\([a-z]*\).png/\1/") $@ 

algscdelta1-2.png algovscdelta1-2.png ovscdelta1-2.png: src1.png src2.png $(RIGID2D)
	$(RIGID2D) src1.png src2.png algscdelta1-2.png scdelta algovscdelta1-2.png ovscdelta1-2.png

algscalar1-2.png algovscalar1-2.png ovscalar1-2.png: src1.png src2.png $(RIGID2D)
	$(RIGID2D) src1.png src2.png algscalar1-2.png scalar algovscalar1-2.png ovscalar1-2.png

algcross1-2.png algovcross1-2.png ovcross1-2.png: src1.png src2.png $(RIGID2D)
	$(RIGID2D) src1.png src2.png algcross1-2.png cross algovcross1-2.png ovcross1-2.png

algscdelta2-1.png algovscdelta2-1.png ovscdelta2-1.png: src1.png src2.png $(RIGID2D)
	$(RIGID2D) src2.png src1.png algscdelta2-1.png scdelta algovscdelta2-1.png ovscdelta2-1.png

algscalar2-1.png algovscalar2-1.png ovscalar2-1.png: src1.png src2.png $(RIGID2D)
	$(RIGID2D) src2.png src1.png algscalar2-1.png scalar algovscalar2-1.png ovscalar2-1.png

algcross2-1.png algovcross2-1.png ovcross2-1.png: src1.png src2.png $(RIGID2D)
	$(RIGID2D) src2.png src1.png algcross2-1.png cross algovcross2-1.png ovcross2-1.png

self-%.png: src1.png src1.png $(GRADMAP)
	$(GRADMAP) src1.png src1.png $(shell echo $@|sed -s "s/self-\([a-z]*\).png/\1/") $@ 

src1grad.png: src1.png $(Image2Grad)
	$(Image2Grad) $< $@

src2grad.png: src2.png $(Image2Grad)
	$(Image2Grad) $< $@

$(RIGID2D) $(Image2Grad) $(GRADMAP): 
	@echo "Error: Please build the project first or set BINROOT in this Makefile to the binary output path"
	@exit 1 
