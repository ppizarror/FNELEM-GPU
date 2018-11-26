/**
FNELEM-GPU - TEST UTILS
Utilitary functions to perform tests.

@package test
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

// Init header file
#ifndef FNELEM_GPU_TEST_FNELEM_UTILS_H
#define FNELEM_GPU_TEST_FNELEM_UTILS_H

// Library imports
#include <cassert>
#include <iostream>
#include <math.h>
#include <string>

/**
 * Check if two numbers are equal under a certain tolerance.
 *
 * @param a Number to evaluate
 * @param b Number to evaluate
 * @return
 */
bool is_num_equal(double a, double b) {
    return fabs(a - b) < 1e-12;
}

/**
 * Print test to console.
 *
 * @param suite Suite name
 * @param test Test name
 */
void test_print_title(const std::string &suite, const std::string &test) {
    std::cout << "[" << suite << "] " << test << std::endl;
}

#endif // FNELEM_GPU_TEST_FNELEM_UTILS_H