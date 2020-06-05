simname = test-sim-6-8-2
fname = $(simname)
SCRATCHDIR = /$(SCRATCH_DIR)/cytof/results/repfam

all: build

build:
	latexmk $(fname).tex -pdf

clean:
	latexmk -c
	rm -f *.bbl *.synctex.gz

cleanall: clean
	latexmk -C

watch:
	latexmk -pvc -pdf $(fname).tex

getresults:
	cd ../../ && \
		make getoutput -s SIM_EXP_NAME=$(simname) && \
		ln -s $(SCRATCHDIR)/$(simname) results
		
