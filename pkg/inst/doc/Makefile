COMPONENTDIR=demos
subsections = \
	$(COMPONENTDIR)/p4.Rnw \
	$(COMPONENTDIR)/p5.Rnw \
	$(COMPONENTDIR)/p15.Rnw \
	$(COMPONENTDIR)/p50.Rnw \
	$(COMPONENTDIR)/p51.Rnw \
	$(COMPONENTDIR)/p55.Rnw \
	$(COMPONENTDIR)/p112.Rnw \
	$(COMPONENTDIR)/p167.Rnw \
	$(COMPONENTDIR)/p376.Rnw \
	$(COMPONENTDIR)/p410.Rnw \
	$(COMPONENTDIR)/p448.Rnw \
	$(COMPONENTDIR)/p489.Rnw \
	$(COMPONENTDIR)/p582.Rnw \
	$(COMPONENTDIR)/p647.Rnw \
	$(COMPONENTDIR)/p660.Rnw \
	$(COMPONENTDIR)/p697.Rnw

all : Companion.pdf ;

Companion.pdf : Companion.tex
	texi2dvi --pdf Companion.tex
  
Companion.tex: $(subsections) Companion.Rnw
	R CMD Sweave Companion.Rnw

%.Rnw : %.H
	cat demos/Sweave_prefix.Rnw $< > $(*F).Rnw
	R CMD Sweave $(*F).Rnw
	cat demos/Stangle_prefix.Rnw $< > $(*F).Rnw
	R CMD Stangle $(*F).Rnw
	grep -v "^###" $(*F).R > ../../demo/$(*F).R
	mv $(*F)*.* $(COMPONENTDIR)
