/**
FNELEM-GPU CPU MATRIX INVERSION
Performs matrix inversion using Gauss Jordan algorithm.

@package fnelem.math
@author ppizarror
@date 20/11/2018
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

// Init header file
#ifndef __FNELEM_MATH_MATRIX_INVERSION_CPU_H
#define __FNELEM_MATH_MATRIX_INVERSION_CPU_H

// Include headers
#include "fematrix.h"

// Constant definition
#define __FEMATRIX_MIN_INVERSION_VALUE 0.0005

/**
 * Performs CPU matrix inversion using Gauss-Jordan elimination algorithm.
 *
 * @param matrix Matrix to inverse
 * @return Inverse matrix
 */
FEMatrix *matrix_inverse_cpu(FEMatrix *matrix);

#endif // __FNELEM_MATH_MATRIX_INVERSION_CPU_H