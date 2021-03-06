SHELL:=/bin/bash

.PHONY: clean getoutput getalloutput killsim metrics postproc run sendoutput
.PHONY: zip_results postproc-and-metrics testrun

# Simulation experiment name
# SIM_EXP_NAME = test-sim-6-3
# SIM_EXP_NAME = test-sim-6-5
# SIM_EXP_NAME = test-sim-6-5-1
# SIM_EXP_NAME = test-sim-6-5-2
# SIM_EXP_NAME = test-sim-6-6
# SIM_EXP_NAME = test-sim-6-6-1
# SIM_EXP_NAME = test-sim-6-7
# SIM_EXP_NAME = test-sim-6-7-1
# SIM_EXP_NAME = test-sim-6-7-2
# SIM_EXP_NAME = test-sim-6-7-3
# SIM_EXP_NAME = test-sim-6-7-4
# SIM_EXP_NAME = test-sim-6-7-5
# SIM_EXP_NAME = test-sim-6-7-6
# SIM_EXP_NAME = test-sim-6-7-7
# SIM_EXP_NAME = test-sim-6-7-8
# SIM_EXP_NAME = test-sim-6-7-9
# SIM_EXP_NAME = test-sim-6-7-10
# SIM_EXP_NAME = test-sim-6-7-12
# SIM_EXP_NAME = test-sim-6-7-13
# SIM_EXP_NAME = test-sim-6-7-14
# SIM_EXP_NAME = test-sim-6-7-15
# SIM_EXP_NAME = test-sim-6-7-16
# SIM_EXP_NAME = test-sim-6-7-17
# SIM_EXP_NAME = test-sim-6-7-18
# SIM_EXP_NAME = test-sim-6-7-19
# SIM_EXP_NAME = test-sim-6-7-20
# SIM_EXP_NAME = test-sim-6-7-21
# SIM_EXP_NAME = test-sim-6-7-22

# NOTE: completed, but unprocessed
# SIM_EXP_NAME = test-sim-6-7-11

SIM_EXP_NAME = test-sim-6-8-1

# Run script
# RUN_SCRIPT = run.jl
# RUN_SCRIPT = run-pt.jl
RUN_SCRIPT = run-imcmc-pt.jl

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
	julia $(RUN_SCRIPT) $(SIM_CONFIG_DIR) $(SIM_RESULTS_DIR) $(SIM_AWS_BUCKET) 0 \
		&> $(SIM_RESULTS_DIR)/master-log.txt &

testrun:
	@echo "Printing results to results/$(SIM_EXP_NAME)-master-log.txt"
	julia $(RUN_SCRIPT) $(SIM_CONFIG_DIR) $(SIM_RESULTS_DIR) $(SIM_AWS_BUCKET) 1 \
		&> results/$(SIM_EXP_NAME)-master-log.txt &

clean:
	rm -rf $(SIM_RESULTS_DIR)/*

getoutput:
	aws s3 sync $(SIM_AWS_BUCKET) $(SIM_RESULTS_DIR) \
		--exclude "*/output.bson" --exclude "*/simdat.bson" \
		--exclude "master-log.txt"

getalloutput:
	aws s3 sync $(SIM_AWS_BUCKET) $(SIM_RESULTS_DIR)

postproc-and-metrics: postproc metrics

postproc:
	julia postprocess.jl $(SIM_RESULTS_DIR)
	make sendoutput -s

sendoutput:
	aws s3 sync $(SIM_RESULTS_DIR) $(SIM_AWS_BUCKET) --exclude '*.nfs*' --exclude '*.zip*'

metrics:
	julia metrics.jl $(SIM_RESULTS_DIR)
	make sendoutput -s

killsim:
	pids=`grep -r pid $(SIM_RESULTS_DIR) | grep -oP '(?<=pid).*' | grep -oP '\d+'`; \
			 for p in $$pids; do kill $$p; done

zip_results:
	cd $(SIM_RESULTS_DIR)/.. && zip -r $(SIM_EXP_NAME).zip $(SIM_EXP_NAME)
	mv $(SIM_RESULTS_DIR).zip results/$(SIM_EXP_NAME).zip

rm_zipped_results:
	rm -i results/$(SIM_EXP_NAME).zip
