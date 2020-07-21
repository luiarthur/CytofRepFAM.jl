#!/bin/bash

jid1=`sbatch submit-job-1.sh`
sbatch --dependency=afterok:$jid1 submit-job-2.sh
