simname = test-sim-6-8-1
fname = $(simname)
SCRATCH_DIR = /scratchdata/$(USER)/cytof/results/repfam

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
		ln -s $(SCRATCH_DIR)/$(simname) results
		
