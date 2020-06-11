simname = run4
SCRATCHDIR = $(SCRATCH_DIR)/cytof/results/repfam

all: build

build:
	latexmk $(simname).tex -pdf

clean:
	latexmk -c
	rm -f *.bbl *.synctex.gz

cleanall: clean
	latexmk -C

watch:
	latexmk -pvc -pdf $(simname).tex

getresults:
	cd ../ && make getresults -s
	ln -s $(SCRATCHDIR)/patients-data/$(simname) results
	Rscript R_plots.R
		
