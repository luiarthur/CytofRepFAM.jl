#!/bin/bash

jid1=$(sbatch submit-job-1.sh | grep -oP '\d+'); echo $jid1
jid2=$(sbatch --dependency=afterok:$jid1 submit-job-2.sh | grep -oP '\d+'); echo $jid2
jid3=$(sbatch --dependency=afterok:$jid2 submit-job-3.sh | grep -oP '\d+'); echo $jid3
