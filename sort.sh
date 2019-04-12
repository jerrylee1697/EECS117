#!/bin/bash
#$ -N mergesort-omp
#$ -q class 
#$ -pe openmp 8

# Module load gcc compiler version 4.8.2
module load  gcc/4.8.2

# NOTE: 
# The OMP_NUM_THREADS variable is automaticaly set to the
# core count of our job.   No need to set it again, but here
# it is.
export OMP_NUM_THREADS=8

# Run OpenMP mergesort on 10 million elements
./mergesort-omp 10000000

