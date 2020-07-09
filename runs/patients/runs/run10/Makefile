.PHONY: run getresults getallresults test-run

simname = run10
data_dir = ../../../data/patients/transformed-data

# NOTE: Using three samples
# data_paths = "$(data_dir)/001_d31_clean.csv,$(data_dir)/007_d35_clean.csv,$(data_dir)/010_d35_clean.csv"

# NOTE: Using two samples
data_paths = "$(data_dir)/001_d31_clean.csv,$(data_dir)/007_d35_clean.csv"

results_dir = /scratchdata/$(USER)/cytof/results/repfam/patients-data/$(simname)
AWS_BUCKET = s3://cytof-repfam/patients-data-results/$(simname)
phi = 0 1 10 100

run:
	for phi in $(phi); do \
		rdir=$(results_dir)/phi$$phi; \
		mkdir -p $$rdir; \
		sleep 3; \
		julia run.jl $$phi $(data_paths) $$rdir $(AWS_BUCKET)/phi$$phi 0 \
			&> $$rdir/log.txt & \
	done

test-run:
	for phi in $(phi); do \
		rdir=$(results_dir)/phi$$phi; \
		mkdir -p $$rdir; \
		sleep 3; \
		julia run.jl $$phi $(data_paths) $$rdir $(AWS_BUCKET)/phi$$phi 1 \
			&> $$rdir/log.txt & \
	done


getallresults:
		aws s3 sync $(AWS_BUCKET) $(results_dir)

getresults:
		aws s3 sync $(AWS_BUCKET) $(results_dir) \
			--exclude "*/output.bson" \
			--exclude "*/log.txt"

sendresults:
		aws s3 sync $(results_dir) $(AWS_BUCKET) \
			--exclude "*/output.bson"

sendallresults:
		aws s3 sync $(results_dir) $(AWS_BUCKET) \
			--exclude '*.nfs'
