/**
FNELEM-GPU - TEST MATH
Force cuda as CPU inversion.

@package test.math
@author ppizarror
@date 22/12/2018
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

#include "../../fnelem/math/matrix_inversion_cpu.h"

/**
 * Performs CPU matrix inversion using Gauss-Jordan elimination algorithm with CUDA, when run
 * by CPU forces CPU function. This is intented to used by CPU tests suite.
 *
 * @param matrix Matrix to inverse
 * @return Inverse matrix
 */
FEMatrix *matrix_inverse_cuda(FEMatrix *feMatrix) {
    return matrix_inverse_cpu(feMatrix);
}