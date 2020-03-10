fname = test-sim-6-7-3

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
		make getoutput -s && \
		make zip_results -s && \
		cd results && unzip $(fname)
