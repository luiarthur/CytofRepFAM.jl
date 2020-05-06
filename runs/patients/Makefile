data_dir = ../data
aws_bucket = s3://cytof-repfam/data/

path_to_transformed_patients_data = $(data_dir)/transformed-patients-data
aws_path_to_transformed_patients_data = $(aws_bucket)/transformed-patients-data

send-transformed-patients-data:
	aws s3 sync $(path_to_transformed_patients_data) $(aws_path_to_transformed_patients_data)

get-transformed-patients-data:
	aws s3 sync $(aws_path_to_transformed_patients_data) $(path_to_transformed_patients_data) 
