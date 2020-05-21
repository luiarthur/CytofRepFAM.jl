.PHONY: run getresults getallresults

data_dir = ../../../data/patients/transformed-data
data_paths = "$(data_dir)/001_d31_clean.csv,$(data_dir)/007_d35_clean.csv,$(data_dir)/010_d35_clean.csv"
results_dir = /scratchdata/$(USER)/cytof/results/repfam/patients-data/run1
AWS_BUCKET = s3://cytof-repfam/patients-data-results/run1

run:
	for phi in "0 1 10 25"; do \
		julia run.jl $$phi $(data_dir) $(results_dir) &; \
		done; \
		aws s3 sync $(results_dir) $(AWS_BUCKET) --exclude '*.nfs'

getallresults:
		aws s3 sync $(AWS_BUCKET) $(results_dir)

getresults:
		aws s3 sync $(AWS_BUCKET) $(results_dir) \
			--exclude "*/output.bson"
