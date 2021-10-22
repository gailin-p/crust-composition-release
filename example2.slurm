#!/bin/bash
  
#SBATCH -J mpi_job
#SBATCH -o mpi_%j.txt
#SBATCH -e mpi_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:10:00

cd /dartfs-hpc/rc/lab/K/Kellerlab/crustal-inversion/runPerplexBatchVp
srun julia mpi_ex.jl
