SHELL:=/bin/bash

.PHONY: cleansim sim retrieveallsim retrievesim postprocsim metrics killsim


# Simulation experiment name
SIM_EXP_NAME = test-sim-6-3

# Simulation settings dir
SIM_CONFIG_DIR = configs/$(SIM_EXP_NAME)

# Scratch directory
SCRATCH_DIR = /scratchdata/$(USER)/cytof/results/repfam

# Simulation results dir
SIM_RESULTS_DIR = $(SCRATCH_DIR)/$(SIM_EXP_NAME)

# Sim AWS bukcet
SIM_AWS_BUCKET = s3://cytof-repfam/$(SIM_EXP_NAME)


### Sim Make cmds ###
sim:
	@echo "Printing results to results/$(SIM_EXP_NAME)-master-log.txt"
	julia run.jl $(SIM_CONFIG_DIR) $(SIM_RESULTS_DIR) $(SIM_AWS_BUCKET) \
		&> results/$(SIM_EXP_NAME)-master-log.txt &

cleansim:
	rm -rf $(SIM_RESULTS_DIR)/*

retrievesim:
	aws s3 sync $(SIM_AWS_BUCKET) $(SIM_RESULTS_DIR) \
		--exclude "*/output.bson" --exclude "*/simdat.bson"

retrieveallsim:
	aws s3 sync $(SIM_AWS_BUCKET) $(SIM_RESULTS_DIR)

postprocsim:
	julia post_process.jl $(SIM_RESULTS_DIR) $(SIM_AWS_BUCKET)

sendsim:
	aws s3 sync $(SIM_RESULTS_DIR) $(SIM_AWS_BUCKET) --exclude '*.nfs*'

metrics:
	julia metrics.jl

notebook:
	jupyter notebook

killsim:
	pids=`grep -r pid $(SIM_RESULTS_DIR) | grep -oP '(?<=pid).*' | grep -oP '\d+'`; \
			 for p in $$pids; do kill $$p; done
