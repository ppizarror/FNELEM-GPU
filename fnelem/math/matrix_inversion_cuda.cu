/**
FNELEM-GPU GPU MATRIX INVERSION
Performs matrix inversion using Gauss Jordan algorithm.
Based on: https://github.com/ZhengzhongSun/Matrix-Inversion-with-CUDA

@package fnelem.math
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
#include <stdio.h>
#include <iostream>
#include "FEMatrix.h"

// Constants
const int MATRIX_INVERSION_CUDA_BLOCKSIZE = 8;

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
 * Save matrix to file.
 *
 * @param A Matrix
 * @param s File name
 * @param n Number of rows
 * @param h Number of columns
 */
void save_matrix_to_file(double *A, std::string s, int n, int h) {
    std::ofstream plik;
    plik.open(s);
    for (int j = 0; j < h; j++) {
        for (int i = 0; i < h; i++) {
            plik << A[j * n + i] << "\t";
        }
        plik << std::endl;
    }
    plik.close();
}

/**
 * Matrix inversion, uses CUDA.
 *
 * @param feMatrix Matrix to inverse
 * @return Inverse matrix
 */
FEMatrix *matrix_inverse_cuda(FEMatrix *feMatrix) {

    // Get matrix
    double *matrix = feMatrix->get_array();

    // Get matrix dimension
    int *matDim = feMatrix->get_dimension();
    int n;
    if (matDim[0] == matDim[1]) {
        n = matDim[0];
    } else {
        throw std::logic_error("Matrix to inverse is not square");
    }

    // Inverse matrix CPU
    double *iMatrix = new double[n * n];

    // Create auxiliar matrices
    double *d_A, *I, *dI;

    // Time of computation
    float time;

    // Create CUDA error handlers
    cudaError_t err;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // Matrix memory size
    int ddsize = n * n * sizeof(double);

    // Creates blocks
    dim3 threadsPerBlock(MATRIX_INVERSION_CUDA_BLOCKSIZE, MATRIX_INVERSION_CUDA_BLOCKSIZE);
    dim3 numBlocks((n + MATRIX_INVERSION_CUDA_BLOCKSIZE - 1) / MATRIX_INVERSION_CUDA_BLOCKSIZE,
                   (n + MATRIX_INVERSION_CUDA_BLOCKSIZE - 1) / MATRIX_INVERSION_CUDA_BLOCKSIZE);

    // Memory allocation
    err = cudaMalloc((void **) &d_A, ddsize);
    if (err != cudaSuccess) {
        std::cout << cudaGetErrorString(err) << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
    }
    err = cudaMalloc((void **) &dI, ddsize);
    if (err != cudaSuccess) {
        std::cout << cudaGetErrorString(err) << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
    }

    // Creates identify matrix
    I = new double[n * n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) I[i * n + i] = 1.0;
            else I[i * n + j] = 0.0;
        }
    }

    // Copy data from CPU to GPU
    err = cudaMemcpy(d_A, matrix, ddsize, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        std::cout << cudaGetErrorString(err) << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
    }
    err = cudaMemcpy(dI, I, ddsize, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        std::cout << cudaGetErrorString(err) << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
    }

    // Timer start
    cudaEventRecord(start, 0);

    // L^(-1)
    for (int i = 0; i < n; i++) {
        nodiag_normalize << < numBlocks, threadsPerBlock >> > (d_A, dI, n, i);
        diag_normalize << < numBlocks, threadsPerBlock >> > (d_A, dI, n, i);
        gaussjordan << < numBlocks, threadsPerBlock >> > (d_A, dI, n, i);
        set_zero << < numBlocks, threadsPerBlock >> > (d_A, dI, n, i);
    }

    // Record cuda events
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    // Copy data from GPU to CPU
    err = cudaMemcpy(iMatrix, dI, ddsize, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        std::cout << cudaGetErrorString(err) << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
    }
    err = cudaMemcpy(I, d_A, ddsize, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        std::cout << cudaGetErrorString(err) << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
    }
    std::cout << "[CUDA] Matrix inversion time: " << time << "ms\n" << std::endl;

    // Free memory
    cudaFree(d_A);
    cudaFree(dI);
    delete[] I;

    // Generate matrix
    return new FEMatrix(iMatrix, n, n);

}