# NOTE: I manu8ally did the following:
# - download data from Google Drivbe
# - Ignore patient all patient 5 data because it is in Excel 2007 format,
#   which I cannot read, as I am on Ubuntu.
# - Due to the file naming convention which uses days, instead of year:
# 	- Changed file name '006_1.5 years_clean.*' to '006_d547_clean.*'
# 	- Changed file name '004_1year_clean.*' to '004_d365_clean.*'
# 	- Changed file name '004_5month_clean.*' to '004_d150_clean.*'
convert:
	python3 sanitize.py
	libreoffice --headless --convert-to csv --outdir data/cutoff data/*.xlsx

