SIM_EXP_NAME = test-sim-6-8-2
awsbucket = s3://cytof-repfam/$(SIM_EXP_NAME)/viz

sendimg:
	aws s3 sync viz $(awsbucket)

getimg:
	aws s3 sync $(awsbucket) viz

# NOTES:
# 1. julia compute_tsne.jl
# 2. Rscript generate_competitor_clusterings.R
# 3. julia get_rgam_best_clustering.jl
# 4. python3 viz.py
