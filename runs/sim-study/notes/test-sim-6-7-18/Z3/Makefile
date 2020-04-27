simname = test-sim-6-7-18
fname = $(simname)-Z3

all: build

build:
	latexmk $(fname).tex -pdf

clean:
	latexmk -c

cleanall: clean
	latexmk -C

watch:
	latexmk -pvc -pdf $(fname).tex

getresults:
	cd ../../../ && \
		make getoutput -s SIM_EXP_NAME=$(simname) && \
		make zip_results -s SIM_EXP_NAME=$(simname) && \
		cd results && unzip $(simname)
