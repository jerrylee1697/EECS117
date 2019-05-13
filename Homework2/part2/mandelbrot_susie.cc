/**
 *  \file mandelbrot_susie.cc
 *
 *  \brief Implement your parallel mandelbrot set in this file.
 */

#include <iostream>
#include <cstdlib>
#include <mpi.h>

int
main (int argc, char* argv[])
{
  /* Lucky you, you get to write MPI code */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,	/* always use this */
			&myrank);	/* process rank, 0 thru N-1 */
    MPI_Comm_size(MPI_COMM_WORLD,	/* always use this */
			&ntasks);	/* #processes in application */
    
    //MPI_SCATTER(sendbuf, sendcount, sendtype, recvbuf, recvcount, 
                      recvtype, root, comm)
 
    // process p
    // P: processes
    // row to compute = p + nP for n = 0, 1, 2
    int n = 0;
    int row = myrank + n * ntasks;
    while (row < height) {
        // Does entire row of p + nP
        y = minY + it * row;
        x = minX;
        for (int j = 0; j < width; ++j) {

            mandelbrot(x, y);
            x += jt;
        }
        n++;
        row = myrank + n * ntasks;
    }
    // https://www.mpi-forum.org/docs/mpi-1.1/mpi-11-html/node70.html
    // https://www.mcs.anl.gov/research/projects/mpi/mpi-standard/mpi-report-1.1/node70.htm
    int *recvbuf;

    if (myrank == 0) {
        MPI_Comm_size(MPI_COMM_WORLD,	/* always use this */
			&ntasks);
            recvbuf = (int *)malloc(width * height * sizeof(int));
    }
    MPI_GATHER(
        sendbuf,
        sendcount,
        sendtype,
        recvbuf,
        recvcount,
        MPI_INT,     /* recvtype */
        0,              /* root */
        MPI_COMM_WORLD  /* comm */
    )

    if (myrank == 0) {
        /*
 * Render all values
 */
        for (i = 0; i < height; ++i) {
            for (j = 0; j < width; ++j) {
                img_view(j, i) = render(results[j][i]/512.0);
            }
        }
        gil::png_write_view("mandelbrot_ms.png", const_view(img));
    }


    MPI_Finalize();		
}

/* eof */
