#!/bin/bash

# Run this script with: sbatch submit-job.sh

#SBATCH -p 128x24   # Partition name
#SBATCH -J cytof-test-run  # Job name
#SBATCH --mail-user=alui2@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH -o out/slurm-job.out    # Name of stdout output file
#SBATCH -N 1         # Total number of nodes requested (128x24/Instructional only)
#SBATCH -n 21        # Total number of mpi tasks requested per node
#SBATCH -t 720:00:00  # Run Time (hh:mm:ss) - 720 hours (optional)
#SBATCH --mem=42G # Memory to be allocated PER NODE

echo "SCRATCH_DIR: $SCRATCH_DIR"

simname="run1"
data_dir="../../../data/patients/transformed-data"
data_paths="${data_dir}/001_d31_clean.csv,${data_dir}/007_d35_clean.csv,${data_dir}/010_d35_clean.csv"
results_dir="${SCRATCH_DIR}/cytof/results/repfam/patients-data/${simname}"
AWS_BUCKET="s3://cytof-repfam/patients-data-results/${simname}"
phi="0 1 10 25"
istest=1

# Load these modules
module load R/R-3.6.1

# Make sure Mclust is installed
echo "Install mclust in R if needed ..."
R -e "install.packages('mclust')"
julia -e 'import Pkg; Pkg.activate("../../../../"); Pkg.build("RCall")'

echo "Doing test run"
for phi in ${phi}
do
  rdir=${results_dir}/phi${phi}
  mkdir -p ${rdir}
	sleep 3
	julia run.jl $phi $data_paths $rdir $AWS_BUCKET/phi$phi $istest
		&> $rdir/log.txt &
done

echo "Done submitting jobs."
