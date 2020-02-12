fname = test-sim-6-5-2

build:
	latexmk $(fname).tex -pdf

clean:
	latexmk -c

cleanall: clean
	latexmk -C

watch:
	latexmk -pvc -pdf $(fname).tex
