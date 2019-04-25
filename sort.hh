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


void mergeSort (keytype* A, int l, int r, int N, keytype* B);

void mergeSort_Serial (keytype* A, int l, int r);

void merge(keytype* A, int l, int m, int r, int N, keytype* B);

void pmergeSort(keytype* A, int p, int r, keytype* B, int s);

void pmerge(keytype* T, int p1, int r1, int p2, int r2, keytype *A, int p3);

void merge_p(keytype* A_start, keytype* A_end, keytype* B_start, keytype* B_end, keytype* R);

int binarySearch(int x, keytype *A, int p, int r);

#endif

/* eof */
