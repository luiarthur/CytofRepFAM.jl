#!/bin/bash

jid1=`sbatch submit-job-1.sh`
sleep 2
sbatch --dependency=afterok:$jid1 submit-job-2.sh
