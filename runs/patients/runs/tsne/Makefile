SHELL:=/bin/bash

imgdir = img
awsbucket = s3://cytof-repfam/patients-data-results/tsne/patients/

all: flowsom mclust tsne sendimg

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

sendimg:
	aws s3 sync img/tsne $(awsbucket)

getimg:
	aws s3 sync $(awsbucket) img/tsne

