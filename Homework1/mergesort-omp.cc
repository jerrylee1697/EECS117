/**
 *  \file mergesort.cc
 *
 *  \brief Implement your mergesort in this file.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <string.h>
// #include <strings.h>

#include "sort.hh"
#include <omp.h>

// using namespace std;

void mySort (int N, keytype* A) {
    keytype* B = newKeys (N);

    #pragma omp parallel shared(A)
    #pragma omp single 
    pmergeSort(A, 0, N-1, B, 0, N);

    memcpy (A, B, N * sizeof (keytype));
    
    free(B);
}

void myParallelMergeSort (int N, keytype* A) {
    #pragma omp parallel 
    #pragma omp single
    mergeSort(A, 0, N-1, N, A);
}

void mySequentialSort (int N, keytype* A) {
    mergeSort_Serial(A, 0, N-1, N, A);
}



void pmergeSort(keytype* A, int p, int r, keytype* B, int s, int size) {
    int n = r - p + 1;
    
    keytype* T = newKeys (n);
    int q = (p + r) / 2;
    int q_1 = q - p + 1;

    if (n == 1) {
        B[s] = A[p];
    }
    else if (n < size/omp_get_max_threads()) { //8192 omp_get_num_procs()
        
        pmergeSort(A, p,     q, T, 0,   size);
        pmergeSort(A, q + 1, r, T, q_1, size);
        pmerge(T, 0, q_1-1, q_1, n-1, B, s);
    }
    else {
        #pragma omp task
        pmergeSort(A, p,     q, T, 0,   size);
        #pragma omp task
        pmergeSort(A, q + 1, r, T, q_1, size);
        #pragma omp taskwait
        pmerge(T, 0, q_1-1, q_1, n-1, B, s);
    }
    free (T);
}


void pmerge(keytype* T, int p1, int r1, int p2, int r2, keytype *A, int p3) {
    int n1 = r1 - p1 + 1;
    int n2 = r2 - p2 + 1;

    if (n1 < n2) {  // Ensure that n1 >= n2
        std::swap(p1, p2);
        std::swap(r1, r2);
        std::swap(n1, n2);
    }
    if (n1 == 0) 
        return;
    if (n1 + n2 <= 8192) {
        merge_p(&T[p1], &T[p1 + n1], &T[p2], &T[p2 + n2], &A[p3]);
    }
    else {
        int q1 = (p1 + r1) / 2;
        int q2 = binarySearch(T[q1], T, p2, r2);
        int q3 = p3 + (q1 - p1) + (q2 - p2);
        A[q3] = T[q1];
        #pragma omp task
        pmerge(T, p1,     q1-1, p2, q2 - 1, A, p3     );
        #pragma omp task
        pmerge(T, q1 + 1, r1,   q2, r2,     A, q3 + 1);
        #pragma omp taskwait
    }
}

void merge_p(keytype* A_start, keytype* A_end, keytype* B_start, keytype* B_end, keytype* R) {
    int counter_a = 0;
    int counter_b = 0;
    int counter_r = 0;
    while (A_start + counter_a < A_end && B_start + counter_b < B_end) {
        if (*(A_start + counter_a) <= *(B_start + counter_b)) {
            *(R + counter_r++) = *(A_start +counter_a++);
        }
        else {
            *(R + counter_r++) = *(B_start + counter_b++);
        }
    }
    while (A_start + counter_a < A_end) {
        *(R + counter_r++) = *(A_start + counter_a++);
    }
    while (B_start  + counter_b < B_end) {
        *(R + counter_r++) = *(B_start + counter_b++);
    }
}


 
void mergeSort (keytype* A, int l, int r, int N, keytype* B) {
    int m = (l + r) / 2;
    if (r-l+1 < N/omp_get_max_threads() && l < r) {
        mergeSort(A, l, m, N, B);
        mergeSort(A, m+1, r, N, B);
        merge(A, l, m, r, N, B);
        return;
    }
    if (l < r) {
        #pragma omp task
        mergeSort(A, l, m, N, B);
        #pragma omp task
        mergeSort(A, m+1, r, N, B);
        #pragma omp taskwait
        merge(A, l, m, r, N, B);
    }
}


void merge(keytype* A, int l, int m, int r, int N, keytype* B) {
    int n1 = m - l + 1;
    int n2 = r - m;

    keytype* L = newKeys (n1);
    keytype* R = newKeys (n2);
    
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
    free(L);
    free(R);
}

int binarySearch(int x, keytype *A, int p, int r) {
    int l = p;
    int h = std::max(p, r + 1);
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

void mergeSort_Serial (keytype* A, int l, int r, int N, keytype* B) {
    if (l < r) {
        int m = (l + r) / 2;
        mergeSort(A, l, m, N, B);
        mergeSort(A, m+1, r, N, B);
        merge(A, l, m, r, N, B);
    }
}

/* eof */