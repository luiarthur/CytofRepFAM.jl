fname = test-sim-6-6

getresults:
	cd ../../ && \
		make getoutput -s && \
		make zip_results -s && \
		cd results && unzip $(fname)

build:
	latexmk $(fname).tex -pdf

clean:
	latexmk -c

cleanall: clean
	latexmk -C

watch:
	latexmk -pvc -pdf $(fname).tex
