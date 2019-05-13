/**
 *  \file mandelbrot_susie.cc
 *
 *  \brief Implement your parallel mandelbrot set in this file.
 */

#include <iostream>
#include <cstdlib>
#include <mpi.h>

int maxRows;    //!< Max number of rows per process

int mandelbrot(double x, double y);

int main (int argc, char* argv[]) {
  /* Lucky you, you get to write MPI code */
    double minX = -2.1;
    double maxX = 0.7;
    double minY = -1.25;
    double maxY = 1.25;

    int height, width;
    if (argc == 3) {
        height = atoi (argv[1]);
        width = atoi (argv[2]);
        assert (height > 0 && width > 0);
    } else {
        fprintf (stderr, "usage: %s <height> <width>\n", argv[0]);
        fprintf (stderr, "where <height> and <width> are the dimensions of the image.\n");
        return -1;
    }

    double it = (maxY - minY)/height;
    double jt = (maxX - minX)/width;
    double x, y;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,	/* always use this */
			&myrank);	/* process rank, 0 thru N-1 */
    MPI_Comm_size(MPI_COMM_WORLD,	/* always use this */
			&ntasks);	/* #processes in application */

    maxRows = ceil(height/ntasks); //!< Sets the maximum number of rows per process
    //MPI_SCATTER(sendbuf, sendcount, sendtype, recvbuf, recvcount, 
     
                    //   recvtype, root, comm)
 
    // process p
    // P: processes
    // row to compute = p + nP for n = 0, 1, 2
    int n = 0;
    int row = myrank + n * ntasks;
    int numberOfRows = 0;
    // http://www.netlib.org/utk/papers/mpi-book/node98.html#SECTION00560000000000000000
    int sendbuf[maxRows * width];   //!< Send buffer from every single process
    int sendbufRow = 0; //!< Keeps track of sendbuf indices

    while (row < height) {
        // Does entire row of p + nP
        y = minY + it * row;
        x = minX;
        for (int j = 0; j < width; ++j) {

            sendbuf[j + sendbufRow * width] = mandelbrot(x, y);
            x += jt;
        }
        n++;
        row = myrank + n * ntasks;
        numberOfRows += 1;
        sendbufRow += 1;
    }
    // https://www.mpi-forum.org/docs/mpi-1.1/mpi-11-html/node70.html
    // https://www.mcs.anl.gov/research/projects/mpi/mpi-standard/mpi-report-1.1/node70.htm
    int *recvbuf;

    if (myrank == 0) {
        MPI_Comm_size(MPI_COMM_WORLD,	/* always use this */
			&ntasks);
            recvbuf = (int *)malloc(width * maxRows * ntasks * sizeof(int));  //!< Receive buffer only by parent process
    }

    MPI_GATHER(
        sendbuf,
        maxRows * width, /* sendcount,*/
        MPI_INT,         /* sendtype,*/ 
        recvbuf,
        maxRows * width, /* recvcount,*/
        MPI_INT,         /* recvtype */
        0,               /* root */
        MPI_COMM_WORLD   /* comm */
    );

    if (myrank == 0) {
        /*
 * Render all values
 */
    gil::rgb8_image_t img(height, width);
    auto img_view = gil::view(img);
        int i, j;
        int bufIndex = 0;
        for (i = 0; i < height; ++i) {
            bufIndex = (i * maxRow * width) * ((ntasks-1) * maxRow * width);
            for (j = 0; j < width; ++j) {
                img_view(j, i) = render(recvbuf(bufIndex + j)/512.0);
            }
        }
        gil::png_write_view("mandelbrot_ms.png", const_view(img));
    }


    MPI_Finalize();		
}

int
mandelbrot(double x, double y) {
  int maxit = 511;
  double cx = x;
  double cy = y;
  double newx, newy;

  int it = 0;
  for (it = 0; it < maxit && (x*x + y*y) < 4; ++it) {
    newx = x*x - y*y + cx;
    newy = 2*x*y + cy;
    x = newx;
    y = newy;
  }
  return it;
}

/* eof */
