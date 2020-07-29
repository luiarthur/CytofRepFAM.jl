simname = run11
fname = $(simname)-phi25
SCRATCHDIR = $(SCRATCH_DIR)/cytof/results/repfam

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
	cd ../ && make getresults -s
	ln -s $(SCRATCHDIR)/patients-data/$(simname) results
	make plots -s

plots:
	cd scripts && Rscript R_plots.R
	cd scripts && python3 more_plots.py	
