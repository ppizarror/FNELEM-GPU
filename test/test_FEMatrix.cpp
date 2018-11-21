/**
FNELEM-GPU - TEST
Test matrix.

@package test
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

// Include sources
#include "../fnelem/math/FEMatrix.h"
#include "../fnelem/math/matrix_inversion_cpu.h"
#include <iostream>
#include <cassert>

void test_matrix_init() {
    FEMatrix matrix = FEMatrix(3, 5);
    matrix.fill_zeros();
    int *dim = matrix.get_dimension();
    assert(dim[0] == 3);
    assert(dim[1] == 5);
    matrix.set(0, 0, 10);
    assert(matrix.get(0, 0) == 10);
}

void test_matrix_disp() {
    FEMatrix matrix = FEMatrix(3, 3);
    matrix.fill_ones();
    matrix.disp();
    assert(matrix.is_square());
}

void test_matrix_cpu_inverse() {
    FEMatrix mat = FEMatrix(3, 3);
    mat.set(0, 0, 1);
    mat.set(0, 1, 2);
    mat.set(0, 2, 3);
    mat.set(1, 0, 5);
    mat.set(1, 1, 2);
    mat.set(1, 2, 1);
    mat.set(2, 0, 2);
    mat.set(2, 1, 2);
    mat.set(2, 2, 3);
    mat.save_to_file("test-matrix-pre-inversion.txt");
    mat.disp();
}

void test_matrix_array() {
    FEMatrix mat = FEMatrix(2, 2);
    double *arr = mat.get_array();
    assert(arr[0] == 0);
}

int test_cpu_inversion() {
    int n = 3;
    auto *L = new double[n * n];
    L[0 * 3 + 0] = 1;
    L[0 * 3 + 1] = 2;
    L[0 * 3 + 2] = 3;
    L[1 * 3 + 0] = 5;
    L[1 * 3 + 1] = 2;
    L[1 * 3 + 2] = 1;
    L[2 * 3 + 0] = 2;
    L[2 * 3 + 1] = 2;
    L[2 * 3 + 2] = 3;
    FEMatrix mat = FEMatrix(L, n, n);
    FEMatrix *matInverse = matrix_inverse_cpu(&mat);
    matInverse->disp();
}

void test_add() {
    FEMatrix m1 = FEMatrix(3, 3);
    m1.fill_ones();
    FEMatrix m2 = FEMatrix(3, 3);
    m2.fill_ones();
    m2.set(0, 0, 3);
    m2.set(1, 1, 5);
    m2.disp();

    // Adds
    m2 += m1;
    assert(m2.get(0, 0) == 4);
    assert(m2.get(1, 1) == 6);

    // Assignation
    m1 = m2;
    assert(m1.get(0, 0) == 4);
    assert(m1.get(1, 1) == 6);
}

void test_substract() {
    FEMatrix m1 = FEMatrix(3, 3);
    m1.fill_ones();
    FEMatrix m2 = FEMatrix(3, 3);
    m2.fill_ones();
    m2.set(0, 0, 3);
    m2.set(1, 1, 5);
    FEMatrix m3 = m2 - m1;
    assert(m3.get(0, 0) == 2);
}

int main() {
    test_matrix_init();
    test_matrix_disp();
    test_matrix_cpu_inverse();
    test_matrix_array();
    test_cpu_inversion();
    test_add();
    test_substract();
    return 0;
}