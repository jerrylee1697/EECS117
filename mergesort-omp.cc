/**
 *  \file mergesort.cc
 *
 *  \brief Implement your mergesort in this file.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "sort.hh"
#include <omp.h>

void
mySort (int N, keytype* A)
{
  /* Lucky you, you get to start from scratch */
    if (N <= 10000) {
        quickSort(N, A);
        return;
    }

    #pragma omp parallel
    #pragma omp single
//   #pragma omp single nowait
    mergeSort(A, 0, N-1, N);
}

 
void mergeSort (keytype* A, int l, int r, int N) {
    if (r-l < N/8) {
        quickSort(r-l+1, (A + l));
        return;
    }
    if (l < r) {
        int m = (l + r) / 2;
        #pragma omp task
        #pragma imp task shared(A)
        mergeSort(A, l, m, N);
        mergeSort(A, m+1, r, N);
        #pragma omp taskwait
        merge(A, l, m, r);
    }
}

void mergeSort_Serial (keytype* A, int l, int r) {
    
    if (l < r) {
        int m = (l + r) / 2;
        mergeSort_Serial(A, l, m);
        mergeSort_Serial(A, m+1, r);
        merge(A, l, m, r);
    }
}

void merge(keytype* A, int l, int m, int r) {
    int n1 = m - l + 1;
    int n2 = r - m;

    keytype* L = new keytype[n1];
    keytype* R = new keytype[n2];
    
    for (int i = 0; i < n1; ++i) {
        L[i] = A[l + i];
    }
    for (int j = 0; j < n2; ++j) {
        R[j] = A[m + 1 + j];
    }

    int i = 0;
    int j = 0; 
    int k = l;

    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            A[k] = L[i];
            ++i;
        }
        else {
            A[k] = R[j];
            ++j;
        }
        ++k;
    }

    while (i <n1) {
        A[k] = L[i];
        ++i;
        ++k;
    }

    while (j <n2) {
        A[k] = R[j];
        ++j;
        ++k;
    }
}

/* eof */