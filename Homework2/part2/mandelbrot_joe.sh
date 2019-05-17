#!/bin/bash
#$ -N Mandelbrot_Joe
#$ -q class
# #$ -pe one-node-mpi 8
#$ -pe mpi 64
#$ -R y

# Grid Engine Notes:
# -----------------
# 1) Use "-R y" to request job reservation otherwise single 1-core jobs
#    may prevent this multicore MPI job from running.   This is called
#    job starvation.

# Module load boost
module load boost/1.57.0

# Module load OpenMPI
# module load openmpi-1.8.3/gcc-4.9.2
module load mpich-3.0.4/gcc-4.8.3

# Run the program 
mpirun -n 8 ./mandelbrot_joe 1000 1000
mpirun -n 16 ./mandelbrot_joe 1000 1000
mpirun -n 32 ./mandelbrot_joe 1000 1000
mpirun -n 64 ./mandelbrot_joe 1000 1000