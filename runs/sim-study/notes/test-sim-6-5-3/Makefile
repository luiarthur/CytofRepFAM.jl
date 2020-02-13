fname = test-sim-6-5-3

build:
	latexmk $(fname).tex -pdf

clean:
	latexmk -c

cleanall: clean
	latexmk -C

watch:
	latexmk -pvc -pdf $(fname).tex
