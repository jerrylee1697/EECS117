/**
 *  \file driver.cc
 *  \brief HW 1: OpenMP
 *
 *  This program
 *
 *  - creates an input array of keys to sort, where the caller gives
 *    the array size as a command-line input;
 *
 *  - sorts it using quicksort, noting the execution time;
 *
 *  - sorts it using YOUR sorting implementation, also noting the
 *    execution time;
 *
 *  - checks that the two sorts produce the same result;
 *
 *  - outputs the execution times and effective sorting rate (i.e.,
 *    keys per second).
 */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "timer.c"
#include <iostream>

#include "sort.hh"
// using namespace std;

/* ============================================================
 */

int
main (int argc, char* argv[])
{
  int N = -1;

  if (argc == 2) {
    N = atoi (argv[1]);
    assert (N > 0);
  } else {
    fprintf (stderr, "usage: %s <n>\n", argv[0]);
    fprintf (stderr, "where <n> is the length of the list to sort.\n");
    return -1;
  }

  stopwatch_init ();
  struct stopwatch_t* timer = stopwatch_create (); assert (timer);

  /* Create an input array of length N, initialized to random values */
  keytype* A_in = newKeys (N);
  for (int i = 0; i < N; ++i)
    A_in[i] = lrand48 ();

  printf ("\nN == %d\n\n", N);

  /* Sort using quicksort */
  keytype* A_qs = newCopy (N, A_in);
  stopwatch_start (timer);
  quickSort (N, A_qs);
  long double t_qs = stopwatch_stop (timer);
  printf ("Quicksort: %Lg seconds ==> %Lg million keys per second\n",
	  t_qs, 1e-6 * N / t_qs);
  assertIsSorted (N, A_qs);


 /* Sort, Sequential Mergesort. */
  keytype* A_s = newCopy (N, A_in);
    keytype* B_s = newCopy (N, A_in);
  stopwatch_start (timer);
  mySequentialSort (N, A_s);
  long double t_s = stopwatch_stop (timer);
  printf ("Sequential Mergesort: %Lg seconds ==> %Lg million keys per second\n",
	  t_s, 1e-6 * N / t_s);
  assertIsSorted (N, A_s);
  assertIsEqual (N, A_s, A_qs);


  /* Sort, Parallel Mergesort. */
  keytype* A_Ps = newCopy (N, A_in);
  keytype* B_Ps = newCopy (N, A_in);
  stopwatch_start (timer);
  myParallelMergeSort (N, A_Ps);
  long double t_Ps = stopwatch_stop (timer);
  printf ("Parallel Mergesort: %Lg seconds ==> %Lg million keys per second\n",
	  t_Ps, 1e-6 * N / t_Ps);
  assertIsSorted (N, A_Ps);
  assertIsEqual (N, A_Ps, A_qs);

 

  /* Sort, calling YOUR routine. */
  keytype* A_ms = newCopy (N, A_in);
  stopwatch_start (timer);
  mySort (N, A_ms);
  long double t_ms = stopwatch_stop (timer);
  printf ("My PMerge sort: %Lg seconds ==> %Lg million keys per second\n",
	  t_ms, 1e-6 * N / t_ms);
  assertIsSorted (N, A_ms);
  assertIsEqual (N, A_ms, A_qs);
  
  /* Cleanup */
  printf ("\n");
  free (A_in);
  free (A_qs);
  free (A_s);
  free (B_s);
  free (A_Ps);
  free (B_Ps);
  free (A_ms);

  stopwatch_destroy (timer);
  return 0;
}

/* eof */
