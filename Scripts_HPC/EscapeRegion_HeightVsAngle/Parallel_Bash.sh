#!/bin/bash

#SBATCH --export=ALL                    # export all environment variables to the batch job
#SBATCH -D .                            # set working directory to .
#SBATCH -p sq                           # submit to the serial queue. Parallel queue is pq
#SBATCH --time=24:00:00                 # maximum walltime for the job
#SBATCH -A Research_Project-T113994     # research project to submit under



#SBATCH -J Testing                  # A single job name for the array
#SBATCH -n 1                        # Number of cores
#SBATCH -N 1                        # All cores on one Node
#SBATCH --mem 10000                # Memory request (10Gb)
#SBATCH -o Testing_%A_%a.out        # Standard output
#SBATCH -e Testing_%A_%a.err        # Standard error

module load Anaconda3

time python Script.py ParamPairs.npz ${SLURM_ARRAY_TASK_ID}
