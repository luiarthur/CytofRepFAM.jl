SHELL:=/bin/bash

.PHONY: clean getoutput getalloutput killsim metrics postproc run sendoutput
.PHONY: zip_results

# Simulation experiment name
# SIM_EXP_NAME = test-sim-6-3
# SIM_EXP_NAME = test-sim-6-5
SIM_EXP_NAME = test-sim-6-5-1

# Simulation settings dir
SIM_CONFIG_DIR = configs/$(SIM_EXP_NAME)

# Scratch directory
SCRATCH_DIR = /scratchdata/$(USER)/cytof/results/repfam

# Simulation results dir
SIM_RESULTS_DIR = $(SCRATCH_DIR)/$(SIM_EXP_NAME)

# Sim AWS bukcet
SIM_AWS_BUCKET = s3://cytof-repfam/$(SIM_EXP_NAME)


### Sim Make cmds ###
run:
	@echo "Printing results to results/$(SIM_EXP_NAME)-master-log.txt"
	julia run.jl $(SIM_CONFIG_DIR) $(SIM_RESULTS_DIR) $(SIM_AWS_BUCKET) \
		&> results/$(SIM_EXP_NAME)-master-log.txt &

clean:
	rm -rf $(SIM_RESULTS_DIR)/*

getoutput:
	aws s3 sync $(SIM_AWS_BUCKET) $(SIM_RESULTS_DIR) \
		--exclude "*/output.bson" --exclude "*/simdat.bson"

getalloutput:
	aws s3 sync $(SIM_AWS_BUCKET) $(SIM_RESULTS_DIR)

postproc:
	@echo "Do post processing ..."
	julia postprocess.jl $(SIM_RESULTS_DIR)
	@echo "Compute metrics ..."
	julia metrics.jl $(SIM_RESULTS_DIR)
	aws s3 sync $(SIM_RESULTS_DIR) $(SIM_AWS_BUCKET) --exclude '*.nfs*'

sendoutput:
	aws s3 sync $(SIM_RESULTS_DIR) $(SIM_AWS_BUCKET) --exclude '*.nfs*'

metrics:
	julia metrics.jl $(SIM_RESULTS_DIR)

killsim:
	pids=`grep -r pid $(SIM_RESULTS_DIR) | grep -oP '(?<=pid).*' | grep -oP '\d+'`; \
			 for p in $$pids; do kill $$p; done

zip_results:
	zip -r results/$(SIM_EXP_NAME).zip $(SIM_RESULTS_DIR)

