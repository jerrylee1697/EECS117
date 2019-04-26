/**
 *  \file sort.hh
 *
 *  \brief Interface to sorting arrays of keys ('keytype' values).
 */

#if !defined (INC_SORT_HH)
#define INC_SORT_HH /*!< sort.hh already included */

/** 'keytype' is the primitive type for sorting keys */
typedef unsigned long keytype;

/**
 *  Sorts an input array containing N keys, A[0:N-1]. The sorted
 *  output overwrites the input array.
 */
void quickSort (int N, keytype* A);

/**
 *  Sorts an input array containing N keys, A[0:N-1]. The sorted
 *  output overwrites the input array. This is the routine YOU will
 *  implement; see 'mergesort-omp.cc'.
 */
void mySort (int N, keytype* A);

/**
 *  Sorts an input array containing N keys, A[0:N-1]. The sorted
 *  output overwrites the input array. This routine uses a parallel
 *  mergesort with a sequential merge and a cutoff for multi-threading.
 */
void myParallelMergeSort (int N, keytype* A);

/**
 *  Sorts an input array containing N keys, A[0:N-1]. 
 *  Uses standard mergesort with no multi-threading.
 */
void mySequentialSort (int N, keytype* A);

/** Returns a new uninitialized array of length N */
keytype* newKeys (int N);

/** Returns a new copy of A[0:N-1] */
keytype* newCopy (int N, const keytype* A);

/**
 *  Checks whether A[0:N-1] is in fact sorted, and if not, aborts the
 *  program.
 */
void assertIsSorted (int N, const keytype* A);

/**
 *  Checks whether A[0:N-1] == B[0:N-1]. If not, aborts the program.
 */
void assertIsEqual (int N, const keytype* A, const keytype* B);

/**
 *  Parallelized mergesort, sequential merge.
 */
void mergeSort (keytype* A, int l, int r, int N, keytype* B);

/**
 *  Standard mergesort + sequential merge.
 */
void mergeSort_Serial (keytype* A, int l, int r, int N, keytype* B);

/**
 *  Sequential merge.
 */
void merge(keytype* A, int l, int m, int r, int N, keytype* B);


/**
 *  Parallelized mergesort. 
 *  Inputs: Unsorted array A, range from p to r
 *          starting index s for array B to store values, size of array
 *  Output: Sorted array B
 */
void pmergeSort(keytype* A, int p, int r, keytype* B, int s, int size);

/**
 *  Parallelized merge
 *  Inputs: Unsorted array T, range p1 to r1, range p2 to r2, 
 *          starting index p3 for where array A stores sorted values
 *  Output: Sorted in array A
 */
void pmerge(keytype* T, int p1, int r1, int p2, int r2, keytype *A, int p3);


/**
 *  Sequential merge used in the parallel merge under a specific cutoff
 *  Outputs sorted array to R
 */
void merge_p(keytype* A_start, keytype* A_end, keytype* B_start, keytype* B_end, keytype* R);

/**
 *  Standard binary search
 */
int binarySearch(int x, keytype *A, int p, int r);

#endif

/* eof */
