#!/bin/bash

jid1=$(sbatch submit-job-1.sh); echo "AAA $jid1 BBB";
jid2=$(sbatch --dependency=afterok:$jid1 submit-job-2.sh); echo $jid2
jid3=$(sbatch --dependency=afterok:$jid2 submit-job-3.sh); echo $jid3
