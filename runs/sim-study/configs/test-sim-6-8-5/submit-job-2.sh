#!/bin/bash

#SBATCH -p 128x24   # Partition name
#SBATCH -J sim-cytof-2  # Job name
#SBATCH --mail-user=alui2@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH -o out/slurm-job.out  # Name of stdout output file
#SBATCH --ntasks=72
#SBATCH --mem-per-cpu=2GB
#SBATCH -t 120:00:00  # Run Time (hh:mm:ss) - 120 hours (optional)

# Run this script with: sbatch submit-job.sh


echo "SCRATCH_DIR: $SCRATCH_DIR"

simname="test-sim-6-8-5"
results_dir="${SCRATCH_DIR}/cytof/results/repfam/${simname}"
aws_bucket="s3://cytof-repfam/${simname}"
phis="0 1 25 100"
zinds="1 2 3"
pmisses="0.0"
istest=0

# Load these modules
module load R/R-3.6.1
module load parallel

# Make sure Mclust is installed
# You can install mclust if needed by first loading the module 
# as done above, and then in (the loaded version of) R:
#`install.packages('mclust')`.
# julia -e 'import Pkg; Pkg.activate(joinpath(@__DIR__, "../../../../")); Pkg.build("RCall")'

echo "Doing test run"
for pmiss in $pmisses
do
  for phi in $phis
  do
    for zind in $zinds
    do
      rdir=${results_dir}/pmiss$pmiss-phi$phi-zind$zind
      mkdir -p ${rdir}
      sleep 3
      # NOTE: The `-n` option specifies the total number of mpi tasks requested per node.
      echo "Current job: $rdir"
      cmd="julia run.jl $results_dir $aws_bucket $pmiss $phi $zind $istest"
      srun -N 1 -n 1 -c 3 --exclusive $cmd &> $rdir/log.txt &
    done
  done
done

echo "Done submitting jobs."
echo "Job submission time:"
date

echo "Jobs are now running. A message will be printed and emailed when jobs are done."

wait

echo "Jobs are completed."
echo "Job completion time:"
date
