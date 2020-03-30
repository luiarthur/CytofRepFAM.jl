fname = test-sim-6-7-7-test

all: build

build:
	latexmk $(fname)-N2000.tex -pdf

clean:
	latexmk -c

cleanall: clean
	latexmk -C

watch:
	latexmk -pvc -pdf $(fname)-N2000.tex

getresults:
	cd ../../ && \
		make getoutput -s SIM_EXP_NAME=$(fname) && \
		make zip_results -s SIM_EXP_NAME=$(fname) && \
		cd results && unzip $(fname)
