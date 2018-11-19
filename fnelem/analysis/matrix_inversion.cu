/**
FNELEM-GPU MATRIX INVERSION

Performs matrix inversion using Gauss Jordan algorithm.
Based on: https://github.com/ZhengzhongSun/Matrix-Inversion-with-CUDA

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

// Library imports
#include <cuda.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string>
#include <vector>

/**
 * NODIAG normalize diagonal matrix (CUDA).
 *
 * @param A Matrix
 * @param I Matrix
 * @param n Dimension
 * @param i Position
 */
__global__ void nodiag_normalize(double *A, double *I, int n, int i) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x < n && y < n)
        if (x == i && x != y) {
            I[x * n + y] /= A[i * n + i];
            A[x * n + y] /= A[i * n + i];
        }
}

/**
 * DIAG normalize diagonal matrix (CUDA).
 *
 * @param A Matrix
 * @param I Matrix
 * @param n Dimension
 * @param i Position
 */
__global__ void diag_normalize(double *A, double *I, int n, int i) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x < n && y < n)
        if (x == y && x == i) {
            I[x * n + y] /= A[i * n + i];
            A[x * n + y] /= A[i * n + i];
        }
}

/**
 * Performs Gauss Jordan algorithm (CUDA).
 *
 * @param A Matrix
 * @param I Matrix
 * @param n Dimension
 * @param i Position
 */
__global__ void gaussjordan(double *A, double *I, int n, int i) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x < n && y < n) {
        if (x != i) {
            I[x * n + y] -= I[i * n + y] * A[x * n + i];
            if (y != i) {
                A[x * n + y] -= A[i * n + y] * A[x * n + i];
            }
        }
    }
}

/**
 * Set zero on matrix (CUDA).
 *
 * @param A Matrix
 * @param I Matrix
 * @param n Dimension
 * @param i Position
 */
__global__ void set_zero(double *A, double *I, int n, int i) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x < n && y < n) {
        if (x != i) {
            if (y == i) {
                A[x * n + y] = 0;
            }
        }
    }
}

/**
 * Matrix inversion, uses CUDA.
 *
 * @param matrix Matrix to inverse
 * @param n Matrix dimension
 * @return Inverse matrix
 */
double *inverse_matrix(double *matrix, int n) {

    // Creates matrices
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
}