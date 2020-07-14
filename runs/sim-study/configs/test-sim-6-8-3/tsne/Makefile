SIM_EXP_NAME = test-sim-6-8-3
awsbucket = s3://cytof-repfam/$(SIM_EXP_NAME)/viz

sendimg:
	aws s3 sync viz $(awsbucket)

getimg:
	aws s3 sync $(awsbucket) viz

graphs:
	julia compute_tsne.jl
	Rscript generate_competitor_clusterings.R
	julia get_rfam_best_clusterings.jl
	python3 viz.py
