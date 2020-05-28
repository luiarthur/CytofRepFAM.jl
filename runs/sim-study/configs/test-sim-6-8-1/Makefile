SCRATCH_DIR = /scratchdata/$(USER)/cytof/results/repfam
SIM_EXP_NAME = test-sim-6-8-1
tsne_dir = $(SCRATCH_DIR)/$(SIM_EXP_NAME)/tsne
aws_dir = s3://cytof-repfam/$(SIM_EXP_NAME)/tsne

sendoutput:
	aws s3 sync $(tsne_dir) $(aws_dir) \
		--exclude '*.nfs*' --exclude '*.zip*'

getoutput:
	aws s3 sync $(aws_dir) $(tsne_dir)

