SHELL = /bin/bash

.PHONY: all reproduce clean-data recreate-all-data
.PHONY: test

### Make commands ###
all: reproduce

# Reproduce results (figures, files, etc.)
reproduce: recreate-all-data

# Recreate data used in CB analysis
# - CB Data (transformed and pre-processed)
recreate-all-data:
	@cd runs/cb && make getcbdata --no-print-directory
	@cd runs/cb && make preproc --no-print-directory

# Remove data used in paper. (In case something goes wrong.)
clean-data:
	@echo "Removing CB and simulated data in `runs/data/`"
	rm -f runs/data/*

# Run tests 
test:
	julia -e 'import Pkg; Pkg.activate("."); Pkg.test();' --color=yes 
