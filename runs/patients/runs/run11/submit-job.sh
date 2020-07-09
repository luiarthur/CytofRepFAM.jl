#!/bin/bash

# Run this script with: sbatch submit-job.sh

#SBATCH -p 128x24   # Partition name
#SBATCH -J run11-cytof  # Job name
#SBATCH --mail-user=alui2@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH -o out/slurm-job.out    # Name of stdout output file
#SBATCH -N 1         # Total number of nodes requested (128x24/Instructional only)
#SBATCH -n 21        # Total number of mpi tasks requested per node
#SBATCH -t 720:00:00  # Run Time (hh:mm:ss) - 720 hours (optional)
#SBATCH --mem=42G # Memory to be allocated PER NODE

echo "SCRATCH_DIR: $SCRATCH_DIR"

simname="run11"
data_dir="../../../data/patients/transformed-data"

# NOTE: using three samples
# data_paths="${data_dir}/001_d31_clean.csv,${data_dir}/007_d35_clean.csv,${data_dir}/010_d35_clean.csv"

# NOTE: using two samples
data_paths="${data_dir}/001_d31_clean.csv,${data_dir}/007_d35_clean.csv"

results_dir="${SCRATCH_DIR}/cytof/results/repfam/patients-data/${simname}"
AWS_BUCKET="s3://cytof-repfam/patients-data-results/${simname}"
phi="0 1 10 100"
istest=0

# Load these modules
module load R/R-3.6.1

# Make sure Mclust is installed
# You can install mclust if needed by first loading the module 
# as done above, and then in (the loaded version of) R:
#`install.packages('mclust')`.
julia -e 'import Pkg; Pkg.activate(joinpath(@__DIR__, "../../../../")); Pkg.build("RCall")'

echo "This is a healthy sign of life ..."

for phi in ${phi}
do
  rdir=${results_dir}/phi${phi}
  mkdir -p ${rdir}
	sleep 3
	julia run.jl $phi $data_paths $rdir $AWS_BUCKET/phi$phi $istest &> $rdir/log.txt &
done

echo "Done submitting jobs."
echo "Job submission time:"
date

echo "Jobs are now running. A message will be printed and emailed when jobs are done."

wait

echo "Jobs are completed."
echo "Job completion time:"
date
