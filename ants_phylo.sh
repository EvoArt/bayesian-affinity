#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=defq                   
#SBATCH --job-name=ants_phylo       
#SBATCH --mail-type=ALL                     
#SBATCH --output=array_%A-%a.log               
#SBATCH --array=1-8

pwd; hostname; date

echo "Running a program on $SLURM_JOB_NODELIST"

module load Workspace/v1

export JULIA_NUM_THREADS=4
julia ants_genus.jl $SLURM_ARRAY_TASK_ID --threads = 4
julia ants_species.jl $SLURM_ARRAY_TASK_ID --threads = 4
                                                                                                            
       