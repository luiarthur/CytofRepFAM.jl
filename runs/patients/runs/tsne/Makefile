SHELL:=/bin/bash

imgdir = img

all: flowsom mclust tsne

flowsom:
	@echo "flowsom: TODO"

plot-mclust-tsne:
	python3 plot_mclust_tsne.py

mclust: preparedata img/mclust-*.csv 

tsne: img/tsne-*.txt

img/tsne-*.txt:
	python3 compute_tsne.py

img/mclust-*.csv:
	Rscript generate_mclust_clusterings.R

preparedata: img/y*-prepped.csv

img/y*-prepped.csv:
	Rscript prepare-data.R

