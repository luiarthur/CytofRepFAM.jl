data_dir = ../data/patients
aws_bucket = s3://cytof-repfam/data/patients

send-data:
	aws s3 sync $(data_dir) $(aws_bucket)

get-data:
	aws s3 sync $(aws_bucket) $(data_dir)
