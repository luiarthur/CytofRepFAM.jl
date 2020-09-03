.PHONY: link_results

simname = run11
RESULTS_DIR = $(SCRATCH_DIR)/cytof/results/repfam/patients-data/$(simname)
AWS_BUCKET = s3://cytof-repfam/patients-data-results/$(simname)

link_results:
	ln -s $(RESULTS_DIR) results

sendimg:
		aws s3 sync $(RESULTS_DIR)/tsne $(AWS_BUCKET)/tsne

getimg:
		aws s3 sync $(AWS_BUCKET)/tsne $(RESULTS_DIR)/tsne

plot:
	python3 make_tsne.py
	Rscript get_flowsom_clust.R
	Rscript get_mclust_clust.R
	python3 make_plots.py
