#!/bin/bash
#$ -q class
#$ -pe mpi 16
#$ -N pingpong
#$ -R y

date
hostname
echo -e "\n\n"

# Module load OpenMPI
# module load openmpi-1.8.3/gcc-4.9.2
module load mpich-3.0.4/gcc-4.8.3

# Run the ping-pong benchmark
mpirun -n 2 ./pingpong

# Generate netplot.png
gnuplot netplot.gnu

# eof
