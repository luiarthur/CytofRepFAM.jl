SHELL = /bin/bash

.PHONY: all  getcbdata
.PHONY: preproc

### Make variables ###
data_dir = ../data

cb_transformed_data_url = https://raw.githubusercontent.com/luiarthur/cytof-data/master/data/cb/cb_transformed.csv

path_to_cb_data = $(data_dir)/cb_transformed.csv

path_to_reduced_cb_data = $(data_dir)/cb_transformed_reduced.csv

# Download transformed CB data.
# For data info, see: https://github.com/luiarthur/cytof-data
getcbdata: $(path_to_cb_data)
	@echo transformed cb data is in $(path_to_cb_data)

$(path_to_cb_data):
	wget $(cb_transformed_data_url) -O $(path_to_cb_data)

# Preprocess transformed CB data
preproc: $(path_to_reduced_cb_data)

$(path_to_reduced_cb_data): $(path_to_cb_data) preprocess_cb_data.jl
	julia preprocess_cb_data.jl
