#!/bin/bash

#SBATCH --time 00-10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20



env | grep SLURM

module load gcc python

DIR=/users/ibarbier
source $DIR/myfirstvenv/bin/activate


python3 TuringInstability.py  #change here


deactivate
