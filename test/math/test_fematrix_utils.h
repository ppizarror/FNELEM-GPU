/**
FNELEM-GPU - FEMATRIX UTILS TEST
Test FEMatrix utils.

@package test.math
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

void __test_fematrix_utils_identity() {
    test_print_title("FEMATRIX-UTILS", "test_fematrix_utils_identity");
    int n = 5;
    FEMatrix *i = FEMatrix_identity(n);
    FEMatrix *it = i->transpose();
    assert(i->is_diag());
    assert(i->sum() == n);
    assert(it->equals(i));
    delete i;
    delete it;
}

void __test_fematrix_utils_vector() {
    test_print_title("FEMATRIX-UTILS", "test_fematrix_utils_vector");
    FEMatrix *vector = FEMatrix_vector(5);
    vector->fill_ones();
    vector->disp();
    std::cout << vector->to_string(true) << std::endl;
    std::cout << vector->to_string(false) << std::endl;
    std::cout << vector->to_string_line() << std::endl;
    delete vector;
}

/**
 * Performs TEST-FEMATRIX-UTILS tests.
 */
void test_fematrix_utils_suite() {
    __test_fematrix_utils_identity();
    __test_fematrix_utils_vector();
}