subsections = \
	p4.tex \
	p5.tex \
	p15.tex \
	p50.tex \
	p51.tex \
	p55.tex \
	p112.tex \
	p167.tex \
	p376.tex \
	p410.tex \
	p448.tex \
	p489.tex \
	p582.tex \
	p647.tex \
	p660.tex \
	p697.tex

all : Companion.pdf clean ;

Companion.pdf : RcompHam94 $(subsections) Companion.tex
	texi2dvi --pdf Companion.tex

Companion.tex: Companion.Rnw
	R CMD Sweave $<

RcompHam94:
	R CMD  build --no-vignettes ../..
	R CMD  INSTALL  --library=.  RcompHam94_0.1.tar.gz

%.tex : %.H
	cat Sweave_prefix.H $< > $(*F).Rnw
	R CMD Sweave $(*F).Rnw
	cat Stangle_prefix.H $< > $(*F).Rnw
	R CMD Stangle $(*F).Rnw
	grep -v "^###" $(*F).R > ../../demo/$(*F).R
	rm $(*F).Rnw

.PHONY : clean 
clean:
	-rm -r RcompHam94 p*.tex *.tar.gz Rplots.pdf *-*.pdf *p*.eps p*.R Companion.toc Companion.out Companion.log Companion.aux Companion.tex



