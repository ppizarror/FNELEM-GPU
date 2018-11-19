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
#include "fnelem/analysis/matrix_utils.cpp"

int main() {

    // Creates matrix
    const int n = 3;
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
    double *iL = inverse_matrix(L, 3);
    save_matrix_to_file(iL, "inverse.txt", n, n);

    return 0;
}