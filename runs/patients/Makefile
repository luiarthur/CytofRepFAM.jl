data_dir = ../data/patients
aws_bucket = s3://cytof-repfam/data/patients

# Send everything in $(data_dir) to S3.
send-data:
	aws s3 sync $(data_dir) $(aws_bucket)

# To get raw, clean, and transformed patients data.
# These will be stored in the dir which $(data_dir) points to
get-data:
	aws s3 sync $(aws_bucket) $(data_dir)
