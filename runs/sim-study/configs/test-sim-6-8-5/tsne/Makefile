SIM_EXP_NAME = test-sim-6-8-5
awsbucket = s3://cytof-repfam/$(SIM_EXP_NAME)/viz

tex_results_dir = "../tex/results"

sendimg:
	aws s3 sync viz $(awsbucket)

getimg:
	aws s3 sync $(awsbucket) viz

graphs:
	julia compute_tsne.jl
	Rscript generate_competitor_clusterings.R
	julia get_rfam_best_clusterings.jl
	python3 viz.py

get-true-zs:
	cp $(tex_results_dir)/pmiss0.0-phi0-zind1/img/Z_true.pdf viz/img/Z1.pdf
	cp $(tex_results_dir)/pmiss0.0-phi0-zind2/img/Z_true.pdf viz/img/Z2.pdf
	cp $(tex_results_dir)/pmiss0.0-phi0-zind3/img/Z_true.pdf viz/img/Z3.pdf
