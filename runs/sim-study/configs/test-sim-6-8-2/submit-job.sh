#!/bin/bash

# Run this script with: sbatch submit-job.sh

#SBATCH -p 128x24   # Partition name
#SBATCH -J sim-cytof  # Job name
#SBATCH --mail-user=alui2@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH -o out/slurm-job.out  # Name of stdout output file
#SBATCH --ntasks=24
#SBATCH -t 240:00:00  # Run Time (hh:mm:ss) - 240 hours (optional)
#SBATCH --mem-per-cpu=2G  # Memory to be allocated per cpu

echo "SCRATCH_DIR: $SCRATCH_DIR"

simname="test-sim-6-8-2"
results_dir="${SCRATCH_DIR}/cytof/results/repfam/${simname}"
aws_bucket="s3://cytof-repfam/${simname}"
phis="0 1 10 25"
zinds="1 2 3"
pmisses="0.0 0.6"
istest=1

# Load these modules
module load R/R-3.6.1

# Make sure Mclust is installed
# You can install mclust if needed by first loading the module 
# as done above, and then in (the loaded version of) R:
#`install.packages('mclust')`.
# julia -e 'import Pkg; Pkg.activate(joinpath(@__DIR__, "../../../../")); Pkg.build("RCall")'

echo "Doing test run"
for phi in $phi
do
  for zind in $zinds
  do
    for pmiss in $pmisses
    do
      rdir=${results_dir}/pmiss$pmiss-phi$phi-zind$zind
      mkdir -p ${rdir}
      sleep 3
      # NOTE: The `-n` option specifies the total number of mpi tasks requested per node.
      srun -N 1 -n 5 --exclusive \
        julia run.jl $results_dir $aws_bucket $pmiss $phi $zin $istest &> $rdir/log.txt &
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
