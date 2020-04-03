fname = test-sim-6-7-9

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
	cd ../../ && \
		make getoutput -s SIM_EXP_NAME=$(fname) && \
		make zip_results -s SIM_EXP_NAME=$(fname) && \
		cd results && unzip $(fname)
