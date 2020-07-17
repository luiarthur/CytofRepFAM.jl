RESULTS_DIR_PREFIX = $(SCRATCH_DIR)/cytof/results/repfam
SIM_EXP_NAME = test-sim-6-8-3
RESULTS_DIR = $(RESULTS_DIR_PREFIX)/$(SIM_EXP_NAME)
AWS_DIR = s3://cytof-repfam/$(SIM_EXP_NAME)

sendoutput:
	aws s3 sync $(RESULTS_DIR) $(AWS_DIR) \
		--exclude '*.nfs*' --exclude '*.zip*'

getalloutput:
	aws s3 sync $(AWS_DIR) $(RESULTS_DIR)

getoutput:
	aws s3 sync $(AWS_DIR) $(RESULTS_DIR) \
		--exclude "*/output.bson" --exclude "*/simdat.bson"  -- exclude "*/log.txt"

