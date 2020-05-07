data_dir = data
aws_bucket = s3://cytof-repfam/data/patients

# NOTE: I manu8ally did the following:
# - download data from Google Drivbe
# - Ignore patient all patient 5 data because it is in Excel 2007 format,
#   which I cannot read, as I am on Ubuntu.
# - Due to the file naming convention which uses days, instead of year:
# 	- Changed file name '006_1.5 years_clean.*' to '006_d547_clean.*'
# 	- Changed file name '004_1year_clean.*' to '004_d365_clean.*'
# 	- Changed file name '004_5month_clean.*' to '004_d150_clean.*'
# - These files are now in my UCSC Google Drive:
#   - https://drive.google.com/drive/u/1/folders/11qiKrKcFNLLYkyeT5pgEVqjSBkYgCWYA
#   - `transformed-data/`: transformed data used in analysis
#   - `data/`: cleaned raw data (consistent marker names and file naming)
#   - `img/`: some pdf of the distribution of each marker expression
#   - Put downloaded data in data/. For example:
#       - data/data/
#       - data/transformed-data/
# - From now on, pull and push data using S3!
convert:
	python3 sanitize.py
	libreoffice --headless --convert-to csv --outdir data/cutoff data/*.xlsx

send-data:
	aws s3 sync $(data_dir) $(aws_bucket)

get-data:
	aws s3 sync $(aws_bucket) $(data_dir)
