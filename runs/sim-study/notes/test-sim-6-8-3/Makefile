simname = test-sim-6-8-3
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
	ln -s $(SCRATCHDIR)/$(simname) results
