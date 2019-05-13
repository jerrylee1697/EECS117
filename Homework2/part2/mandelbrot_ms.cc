/**
 *  \file mandelbrot_ms.cc
 *
 *  \brief Implement your parallel mandelbrot set in this file.
 */


/*
 *# Run the program 
    mpirun -np 64  ./mandelbrot_ms 1000 1000
    -np <number of processes>
 * 
 */
#include <iostream>
#include <cstdlib>
#include <assert.h> 
#include "render.hh" 

using namespace std;
// int
// main (int argc, char* argv[])
// {
//   /* Lucky you, you get to write MPI code */
// }

/* eof */

#include <mpi.h>

#define WORKTAG		1
#define DIETAG		2

void master();
void slave();
int mandelbrot(double x, double y);

int height;
int width;

int main (int argc, char* argv[]) {
    if (argc == 3) {
        height = atoi (argv[1]);
        width = atoi (argv[2]);
        assert (height > 0 && width > 0);
    } else {
        fprintf (stderr, "usage: %s <height> <width>\n", argv[0]);
        fprintf (stderr, "where <height> and <width> are the dimensions of the image.\n");
        return -1;
    }


	int		myrank;

	MPI_Init(&argc, &argv);		/* initialize MPI */
	MPI_Comm_rank(MPI_COMM_WORLD,	/* always use this */
			&myrank);	/* process rank, 0 thru N-1 */

	if (myrank == 0) {
		master();
	} else {
		slave();
	}

	MPI_Finalize();			/* cleanup MPI */
}

void master() {
	int		ntasks, rank, work;
	int		result;
	MPI_Status	status;
    
    double minX = -2.1;
    double maxX = 0.7;
    double minY = -1.25;
    double maxY = 1.25;

    double it = (maxY - minY)/height;
    double jt = (maxX - minX)/width;
    double x, y;


    gil::rgb8_image_t img(height, width);
    auto img_view = gil::view(img);

	MPI_Comm_size(MPI_COMM_WORLD,	/* always use this */
			&ntasks);	/* #processes in application */
/*
 * Seed the slaves.
 */
    // height, width, {x,y} to be passed
    int values[width][height][2];
    int results[width][height];

    y = minY;
    for (int i = 0; i < height; ++i) {
        x = minX;
        for (int j = 0; j < width; ++j) {
            // rank = i * width + j;
            values[j][i][0] = x;
            values[j][i][1] = y;
            results[j][i] = 0;
            // MPI_Send(values[i][j],		/* message buffer */
			// 2,		/* one data item */
			// MPI_INT,	/* data item is an integer */
			// rank,		/* destination process rank */
			// WORKTAG,	/* user chosen message tag */
			// MPI_COMM_WORLD);/* always use this */
            x += jt;
        }
        y += it;
    }

    int i = 0;
    int j = 0;
	for (rank = 1; rank < ntasks; ++rank) {
        i = i + j/width;
        j = j % width;

		// work = /* get_next_work_request */;

		MPI_Send(values[j][i],		/* message buffer */
			2,		/* one data item */
			MPI_INT,	/* data item is an integer */
			rank,		/* destination process rank */
			WORKTAG,	/* user chosen message tag */
			MPI_COMM_WORLD);/* always use this */
        j++;
	}
    std::cout << "Enters while\n";
/*
 * Receive a result from any slave and dispatch a new work request
 * work requests have been exhausted.
 */
	// work = /* get_next_work_request */;
    int i_recv = 0;
    int j_recv = 0;

	while (i * width + j < width * height) {
        i = i + j/width;
        j = j % width;
        
        i_recv = i_recv + j_recv / width;
        j_recv = j_recv % width;

		MPI_Recv(&result,	/* message buffer */
			1,		/* one data item */
			MPI_INT,	/* data item is a double real */
			MPI_ANY_SOURCE,	/* receive from any sender */
			MPI_ANY_TAG,	/* receive any type of message */
			MPI_COMM_WORLD,	/* always use this */
			&status);	/* info about received message */
        results[j_recv][i_recv] = result;
        j_recv++;

		MPI_Send(values[j][i], 2, MPI_INT, status.MPI_SOURCE,
				WORKTAG, MPI_COMM_WORLD);

		// work = /* get_next_work_request */;
        j++;
	}
    std::cout << "Exits while\n";
/*
 * Receive results for outstanding work requests.
 */
	for (rank = 1; rank < ntasks; ++rank) {
        i_recv = i_recv + j_recv/width;
        j_recv = j_recv%width;
		MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE,
			    MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        results[j_recv][i_recv] = result;
        j_recv++;
	}
    std::cout << "Finishes gathering "<< i_recv << ' ' << j_recv << std::endl;;

/*
 * Render all values
 */
    for (i = 0; i < height; ++i) {
        for (j = 0; j < width; ++j) {
            img_view(j, i) = render(results[j][i]/512.0);
        }
    }
    gil::png_write_view("mandelbrot_ms.png", const_view(img));

/*
 * Tell all the slaves to exit.
 */
	for (rank = 1; rank < ntasks; ++rank) {
		MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
	}
}

void slave() {
	int		result;
	int		    work[2];
	MPI_Status	status;

	for (;;) {
		MPI_Recv(work, 2, MPI_INT, 0, MPI_ANY_TAG,
				MPI_COMM_WORLD, &status);
        int x = work[0];
        int y = work[1];
        std::cout << "Working on " << x << ' ' << y << std::endl;

/*
 * Check the tag of the received message.
 */
		if (status.MPI_TAG == DIETAG) {
			return;
		}

		result = mandelbrot(x,y);

		MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
}

int mandelbrot(double x, double y) {
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