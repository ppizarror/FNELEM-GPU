/**
FNELEM-GPU - FEMATRIX UTILS TEST
Test FEMatrix utils.

@package test.fnelem.math
@author ppizarror
@date 21/11/2018
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
#include "../test_utils.h"
#include "../../fnelem/math/fematrix.h"
#include "../../fnelem/math/fematrix_utils.h"

void test_identity() {
    int n = 5;
    FEMatrix i = FEMatrix_identity(n);
    assert(i.is_diag());
    assert(i.sum() == n);
    assert(i.transpose() == i);
}

void test_vector() {
    FEMatrix vector = FEMatrix_vector(5);
    vector.fill_ones();
    vector.disp();
}

int main() {
    test_identity();
    test_vector();
    return 0;
}