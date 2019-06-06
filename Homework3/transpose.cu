#include <stdlib.h>
#include <stdio.h>

#include "cuda_utils.h"
#include "timer.c"

typedef float dtype;


__global__ 
void matTrans(dtype* AT, dtype* A, int N)  {
	/* Fill your code here */
    __shared__ float tile[32][32+1];
    
    int x = blockIdx.x * 32 + threadIdx.x;
    int y = blockIdx.y * 32 + threadIdx.y;
    int width = gridDim.x * 32;
  
    for (int j = 0; j < 32; j += 8)
       tile[threadIdx.y+j][threadIdx.x] = A[(y+j)*width + x];
  
    __syncthreads();
  
    x = blockIdx.y * 32 + threadIdx.x;  // transpose block offset
    y = blockIdx.x * 32 + threadIdx.y;
  
    for (int j = 0; j < 32; j += 8) {
       AT[(y+j)*width + x] = tile[threadIdx.x][threadIdx.y + j];
    }
}

void
parseArg (int argc, char** argv, int* N)
{
	if(argc == 2) {
		*N = atoi (argv[1]);
		assert (*N > 0);
	} else {
		fprintf (stderr, "usage: %s <N>\n", argv[0]);
		exit (EXIT_FAILURE);
	}
}


void
initArr (dtype* in, int N)
{
	int i;

	for(i = 0; i < N; i++) {
		in[i] = (dtype) rand () / RAND_MAX;
	}
}

void
cpuTranspose (dtype* A, dtype* AT, int N)
{
	int i, j;

	for(i = 0; i < N; i++) {
		for(j = 0; j < N; j++) {
			AT[j * N + i] = A[i * N + j];
		}
	}
}

int
cmpArr (dtype* a, dtype* b, int N)
{
	int cnt, i;

	cnt = 0;
	for(i = 0; i < N; i++) {
		if(abs(a[i] - b[i]) > 1e-6) cnt++;
	}

	return cnt;
}



void
gpuTranspose (dtype* A, dtype* AT, int N)
{
    struct stopwatch_t* timer = NULL;
    long double t_gpu;

    int pad = 0;
    if (N%32 != 0) {
        pad = 32 - N % 32;
    }
    dim3 dimGrid((N + pad)/32, (N + pad)/32, 1);
    dim3 dimBlock(32, 8, 1);
    // fprintf (stderr, "Finish dim3s\n");

    // Create temp in for padding
    dtype *tempIn = (dtype*) malloc ((N + pad) * (N + pad) * sizeof (dtype));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            tempIn[i*(N + pad) + j] = A[i * N + j];
        }
    }
    
    /* Cuda malloc*/
    dtype *idata, *tdata;
    cudaMalloc(&idata, (N + pad) * (N + pad) * sizeof (dtype));
    cudaMemcpy(idata, tempIn, (N + pad) * (N + pad) * sizeof (dtype), cudaMemcpyHostToDevice);
    // fprintf (stderr,  "Finish GPU Mallocs\n");
    cudaMalloc(&tdata, (N + pad) * (N + pad) * sizeof (dtype));
    // fprintf (stderr,  "Finish Memcopy1\n");
    /* Setup timers */
    stopwatch_init ();
    timer = stopwatch_create ();

    stopwatch_start (timer);
    /* run your kernel here */
    matTrans<<<dimGrid, dimBlock>>>(tdata, idata, N + pad);
  
    cudaThreadSynchronize ();
    t_gpu = stopwatch_stop (timer);
    // fprintf (stderr,  "Finish matrix Trans\n");
    
    // Undo padding
    dtype* tempOut = (dtype*) malloc ((N + pad) * (N + pad) * sizeof (dtype));
    cudaMemcpy(tempOut, tdata, (N + pad) * (N + pad) * sizeof (dtype), cudaMemcpyDeviceToHost);
    for (int i = 0; i < N * N; ++i) {
        for (int j = 0; j < N * N; ++j) {
            AT[i * N + j] = tempOut[i * (N + pad) + j];
        }
    }
    free (tempOut);
    free (tempIn);


    fprintf (stderr, "Size N: %d \n", N);
    fprintf (stderr, "GPU transpose: %Lg secs ==> %Lg billion elements/second\n",
           t_gpu, (N * N) / t_gpu * 1e-9 );

    cudaFree(idata);
    cudaFree(tdata);
    // fprintf (stderr, "Free\n");
}

int 
main(int argc, char** argv)
{
  /* variables */
	dtype *A, *ATgpu, *ATcpu;
  int err;

	int N;

  struct stopwatch_t* timer = NULL;
  long double t_cpu;


	N = -1;
	parseArg (argc, argv, &N);

  /* input and output matrices on host */
  /* output */
  ATcpu = (dtype*) malloc (N * N * sizeof (dtype));
  ATgpu = (dtype*) malloc (N * N * sizeof (dtype));

  /* input */
  A = (dtype*) malloc (N * N * sizeof (dtype));

	initArr (A, N * N);

	/* GPU transpose kernel */
	gpuTranspose (A, ATgpu, N);

  /* Setup timers */
  stopwatch_init ();
  timer = stopwatch_create ();

	stopwatch_start (timer);
  /* compute reference array */
	cpuTranspose (A, ATcpu, N);
  t_cpu = stopwatch_stop (timer);
  fprintf (stderr, "Time to execute CPU transpose kernel: %Lg secs\n",
           t_cpu);

  /* check correctness */
	err = cmpArr (ATgpu, ATcpu, N * N);
	if(err) {
		fprintf (stderr, "Transpose failed: %d\n", err);
	} else {
		fprintf (stderr, "Transpose successful\n");
	}

	free (A);
	free (ATgpu);
	free (ATcpu);

  return 0;
}
