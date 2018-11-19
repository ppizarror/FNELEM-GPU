/**
FNELEM-GPU MAIN FILE
Performs finite element structural analysis using an 4-node membrane, matrix inversion
was calculated using a CUDA algorithm (Gauss Jordan inversion).

@package fnelem.analysis
@author ppizarror
@date 19/11/2018
@license
    MIT License
    Copyright (c) 2018 Pablo Pizarro R.

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

// CUDA library imports
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>

// FNELEM library imports
#include "fnelem/analysis/matrix_inversion.cu"

int main() {

    const int n = 3;
    // creating input
    double *iL = new double[n * n];
    double *L = new double[n * n];

    L[0 * 3 + 0] = 1;
    L[0 * 3 + 1] = 2;
    L[0 * 3 + 2] = 3;
    L[1 * 3 + 0] = 5;
    L[1 * 3 + 1] = 2;
    L[1 * 3 + 2] = 1;
    L[2 * 3 + 0] = 2;
    L[2 * 3 + 1] = 2;
    L[2 * 3 + 2] = 3;

    // matrix_read(L, n);
    save_inverse_matrix_to_file(L, "inv4.txt", n, n);
    //savetofile(L, "L.txt", n, n);

    cout << "inv\n";
    double *d_A, *d_L, *I, *dI;
    float time;
    cudaError_t err;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    int ddsize = n * n * sizeof(double);
    int blocksize = 8;

    dim3 threadsPerBlock(blocksize, blocksize);
    dim3 numBlocks((n + blocksize - 1) / blocksize, (n + blocksize - 1) / blocksize);
    // memory allocation
    err = cudaMalloc((void **) &d_A, ddsize);
    if (err != cudaSuccess) {
        cout << cudaGetErrorString(err) << " in " << __FILE__ << " at line " << __LINE__ << endl;
    }
    err = cudaMalloc((void **) &dI, ddsize);
    if (err != cudaSuccess) {
        cout << cudaGetErrorString(err) << " in " << __FILE__ << " at line " << __LINE__ << endl;
    }
    I = new double[n * n];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) I[i * n + i] = 1.0;
            else I[i * n + j] = 0.0;
        }
    }

    //copy data from CPU to GPU
    err = cudaMemcpy(d_A, L, ddsize, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        cout << cudaGetErrorString(err) << " in " << __FILE__ << " at line " << __LINE__ << endl;
    }
    err = cudaMemcpy(dI, I, ddsize, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        cout << cudaGetErrorString(err) << " in " << __FILE__ << " at line " << __LINE__ << endl;
    }

    //timer start
    cudaEventRecord(start, 0);

    // L^(-1)
    for (int i = 0; i < n; i++) {
        nodiag_normalize << < numBlocks, threadsPerBlock >> > (d_A, dI, n, i);
        diag_normalize << < numBlocks, threadsPerBlock >> > (d_A, dI, n, i);
        gaussjordan << < numBlocks, threadsPerBlock >> > (d_A, dI, n, i);
        set_zero << < numBlocks, threadsPerBlock >> > (d_A, dI, n, i);
    }

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    //copy data from GPU to CPU
    err = cudaMemcpy(iL, dI, ddsize, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        cout << cudaGetErrorString(err) << " in " << __FILE__ << " at line " << __LINE__ << endl;
    }
    err = cudaMemcpy(I, d_A, ddsize, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        cout << cudaGetErrorString(err) << " in " << __FILE__ << " at line " << __LINE__ << endl;
    }

    cout << "Cuda Time - inverse: " << time << "ms\n";
    save_inverse_matrix_to_file(iL, "inv1.txt", n, n);
    save_inverse_matrix_to_file(I, "inv2.txt", n, n);
    save_inverse_matrix_to_file(L, "inv3.txt", n, n);
    //savetofile(I, "I.txt", n, n);
    //savetofile(I, "I.txt", n, n);
    cudaFree(d_A);
    cudaFree(dI);

    double *c = new double[n * n];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            c[i * n + j] = 0;  //put the initial value to zero
            for (int x = 0; x < n; x++)
                c[i * n + j] = c[i * n + j] + L[i * n + x] * iL[x * n + j];  //matrix multiplication
        }
    save_inverse_matrix_to_file(c, "c.txt", n, n);

    delete[]I;
    delete[]L;
    delete[]iL;













    /**
    const int arraySize = 5;
    const int a[arraySize] = {1, 2, 3, 4, 5};
    const int b[arraySize] = {10, 20, 30, 40, 50};
    int c[arraySize] = {0};

    // Add vectors in parallel.
    cudaError_t cudaStatus = addWithCuda(c, a, b, arraySize);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");
        return 1;
    }

    printf("{1,2,3,4,5} + {10,20,30,40,50} = {%d,%d,%d,%d,%d}\n",
           c[0], c[1], c[2], c[3], c[4]);

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }

    return 0;
}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size) {
    int *dev_a = 0;
    int *dev_b = 0;
    int *dev_c = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void **) &dev_c, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void **) &dev_a, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void **) &dev_b, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    addKernel << < 1, size >> > (dev_c, dev_a, dev_b);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }

    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    Error:
    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);
    */
    return 0;
}