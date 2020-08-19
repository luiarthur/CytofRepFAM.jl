simname = test-sim-6-8-7
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
	rm -f results
	ln -s $(SCRATCHDIR)/$(simname) results
	make plots -s 

plots:
	cd scripts && python3 make_plots.py
