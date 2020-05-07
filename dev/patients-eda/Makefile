data_dir = data
aws_bucket = s3://cytof-repfam/data/patients

send-data:
	aws s3 sync $(data_dir) $(aws_bucket) --exclude "*/.gitignore"

get-data:
	aws s3 sync $(aws_bucket) $(data_dir)

# NOTE: I manually did the following:
# - download data from Google Drive. Put in `data/raw`
# - Ignore patient all patient 5 data because it is in Excel 2007 format,
#   which I cannot read, as I am on Ubuntu.
# - Due to the file naming convention which uses days, instead of year:
# 	- Changed file name '006_1.5 years_clean.*' to '006_d547_clean.*'
# 	- Changed file name '004_1year_clean.*' to '004_d365_clean.*'
# 	- Changed file name '004_5month_clean.*' to '004_d150_clean.*'
# - These files are now in my UCSC Google Drive:
#   - https://drive.google.com/drive/u/1/folders/11qiKrKcFNLLYkyeT5pgEVqjSBkYgCWYA
#   - `data/transformed-data/`: transformed data used in analysis
#   - `data/raw`: cleaned raw data (consistent marker names and file naming)
#   - `data/cutoff`: cleaned cutoff files in csv format
#   - `img/`: some pdf of the distribution of each marker expression
convert:
	python3 sanitize.py
	libreoffice --headless --convert-to csv --outdir data/cutoff data/raw/*.xlsx

zip-data:
	zip -r data.zip data
