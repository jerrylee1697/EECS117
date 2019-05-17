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
double t_start, t_elapsed;      /* Timer variables */

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

    // double t_start, t_elapsed;      /* Timer variables */

	MPI_Init(&argc, &argv);		/* initialize MPI */
	MPI_Comm_rank(MPI_COMM_WORLD,	/* always use this */
			&myrank);	/* process rank, 0 thru N-1 */

    t_start = MPI_Wtime ();         /* Start Timer */
    

	if (myrank == 0) {
		master();
	} else {
		slave();
	}

	MPI_Finalize();			/* cleanup MPI */
}

void master() {
	int		ntasks, rank, work;
	int		result[3];  // Receives value, index j, index i
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
    double values[width][height][4];
    int results[width][height];

    y = minY;
    for (int i = 0; i < height; ++i) {
        x = minX;
        for (int j = 0; j < width; ++j) {
            // rank = i * width + j;
            values[j][i][0] = x;
            values[j][i][1] = y;
            values[j][i][2] = j;    //index j
            values[j][i][3] = i;    //index i
            results[j][i] = 0;
            x += jt;
        }
        y += it;
    }

    int i = 0;  // Value counter Height
    int j = 0;  // Value counter Width
	for (rank = 1; rank < ntasks; ++rank) {
        

		// work = /* get_next_work_request */;

		MPI_Send(values[j][i],		/* message buffer */
			4,		/* one data item */
			MPI_DOUBLE,	/* data item is an integer */
			rank,		/* destination process rank */
			WORKTAG,	/* user chosen message tag */
			MPI_COMM_WORLD);/* always use this */

        j++;

        i = i + j/width;
        j = j % width;
	}
    // std::cout << "Enters while\n";
/*
 * Receive a result from any slave and dispatch a new work request
 * work requests have been exhausted.
 */
	// work = /* get_next_work_request */;
    int i_recv = 0;
    int j_recv = 0;

	// while (i * width + j < width * height) {
    while (i < height && j < width) {

		MPI_Recv(result,	/* message buffer */
			3,		        /* 3 data items: value, index j, index i */
			MPI_INT,	    /* data item is a double real */
			MPI_ANY_SOURCE,	/* receive from any sender */
			MPI_ANY_TAG,	/* receive any type of message */
			MPI_COMM_WORLD,	/* always use this */
			&status);	    /* info about received message */
            // [index j][index i]
        results[result[1]][result[2]] = result[0];

		MPI_Send(values[j][i], 4, MPI_DOUBLE, status.MPI_SOURCE,
				WORKTAG, MPI_COMM_WORLD);

        j++;

        j_recv++;

        i = i + j/width;
        j = j % width;
        
        i_recv = i_recv + j_recv / width;
        j_recv = j_recv % width;
	}
    // std::cout << "Exits while\n";
/*
 * Receive results for outstanding work requests.
 */
	for (rank = 1; rank < ntasks; ++rank) {
        

		MPI_Recv(result, 3, MPI_INT, MPI_ANY_SOURCE,
			    MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        results[result[1]][result[2]] = result[0];

        j_recv++;

        i_recv = i_recv + j_recv / width;
        j_recv = j_recv % width;
	}
    // std::cout << "Finishes gathering "<< i_recv << ' ' << j_recv << std::endl;;

    
/*
 * Render all values
 */
    for (i = 0; i < height; ++i) {
        for (j = 0; j < width; ++j) {
            img_view(j, i) = render(results[j][i]/512.0);
        }
    }
    gil::png_write_view("mandelbrot_ms.png", const_view(img));
    t_elapsed = MPI_Wtime () - t_start;         /* Get end time */
    cout << "ntasks: " << ntasks << " time: "<< t_elapsed << endl;

/*
 * Tell all the slaves to exit.
 */
	for (rank = 1; rank < ntasks; ++rank) {
		MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
	}
}

void slave() {
	int		    result[3];
	double		work[4];
	MPI_Status	status;

	for (;;) {
		MPI_Recv(work, 4, MPI_DOUBLE, 0, MPI_ANY_TAG,
				MPI_COMM_WORLD, &status);
        
        // std::cout << "Working on " << x << ' ' << y << std::endl;

/*
 * Check the tag of the received message.
 */
		if (status.MPI_TAG == DIETAG) {
			return;
		}
        double x = work[0];
        double y = work[1];

		result[0] = mandelbrot(x,y);
        result[1] = (int)work[2];
        result[2] = (int)work[3];

		MPI_Send(result, 3, MPI_INT, 0, 0, MPI_COMM_WORLD);
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