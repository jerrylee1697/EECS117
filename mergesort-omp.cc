/**
 *  \file mergesort.cc
 *
 *  \brief Implement your mergesort in this file.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "sort.hh"
#include <omp.h>

using namespace std;

void
mySort (int N, keytype* A)
{
  /* Lucky you, you get to start from scratch */
    // if (N <= 10000) {
    //     quickSort(N, A);
    //     return;
    // }

    
//   #pragma omp single nowait
    // mergeSort(A, 0, N-1, N);
    keytype* B = newKeys (N);
    #pragma omp parallel
    #pragma omp single
    pmergeSort(A, 0, N-1, B, 0);
    A = B;
}

 
void mergeSort (keytype* A, int l, int r, int N) {
    if (r-l < N/8) {
        quickSort(r-l+1, (A + l));
        return;
    }
    if (l < r) {
        int m = (l + r) / 2;
        #pragma omp task
        #pragma omp task shared(A)
        mergeSort(A, l, m, N);
        mergeSort(A, m+1, r, N);
        #pragma omp taskwait
        merge(A, l, m, r);
    }
}

void pmergeSort(keytype* A, int p, int r, keytype* B, int s) {
    int n = r - p + 1;
    if (n == 1) {
        B[s] = A[p];
    }
    else {
        keytype* T = newKeys (n);
        int q = (p + r) / 2;
        int q_1 = q - p + 1;
        #pragma omp task
        pmergeSort(A, p, q, T, 1);
        pmergeSort(A, q + 1, r, T, q_1 + 1);
        #pragma omp taskwait
        pmerge(T, 1, q_1, q_1 + 1, n, B, s);
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

void pmerge(keytype* T, int p1, int r1, int p2, int r2, keytype *A, int p3) {
    int n1 = r1 - p1 + 1;
    int n2 = r2 - p2 + 1;

    if (n1 < n2) {  // Ensure that n1 >= n2
        swap(p1, p2);
        swap(r1, r2);
        swap(n1, n2);
        // p1 = p2;
        // r1 = r2;
        // n1 = n2;
    }
    if (n1 == 0) 
        return;
    else {
        int q1 = (p1 + r1) / 2;
        int q2 = binarySearch(A[q1], A, p2, r2);
        int q3 = p3 + (q1 - p1) + (q2 - p2);
        A[q3] = T[q1];
        #pragma omp task
        pmerge(T, p1, q1-1, p2, q2 - 1, A, p3);
        pmerge(T, q1 + 1, r1, q2, r2, A, q3 + 1);
        #pragma omp taskwait
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

int binarySearch(int x, keytype *A, int p, int r) {
    int l = p;
    int h = max(p, r + 1);
    while (l < h) {
        int mid = (l + h) / 2;
        if (x <= A[mid]) {
            h = mid;
        }
        else {
            l = mid + 1;
        }
    }
    return h;
}
/* eof */