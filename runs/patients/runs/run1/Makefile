.PHONY: run getresults getallresults test-run

data_dir = ../../../data/patients/transformed-data
data_paths = "$(data_dir)/001_d31_clean.csv,$(data_dir)/007_d35_clean.csv,$(data_dir)/010_d35_clean.csv"
results_dir = /scratchdata/$(USER)/cytof/results/repfam/patients-data/run1
AWS_BUCKET = s3://cytof-repfam/patients-data-results/run1
phi = 0 1 10 25

run:
	for phi in $(phi); do \
		rdir=$(results_dir)/phi$$phi; \
		mkdir -p $$rdir; \
		sleep 3; \
		julia run.jl $$phi $(data_paths) $$rdir $(AWS_BUCKET)/$$rdir 0 \
			&> $$rdir/log.txt & \
	done

test-run:
	for phi in $(phi); do \
		rdir=$(results_dir)/phi$$phi; \
		mkdir -p $$rdir; \
		sleep 3; \
		julia run.jl $$phi $(data_paths) $$rdir $(AWS_BUCKET)/$$rdir 1 \
			&> $$rdir/log.txt & \
	done


getallresults:
		aws s3 sync $(AWS_BUCKET) $(results_dir)

getresults:
		aws s3 sync $(AWS_BUCKET) $(results_dir) \
			--exclude "*/output.bson"
