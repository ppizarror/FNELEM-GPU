/**
FNELEM-GPU - TEST LOAD PATTERN CONSTANT
Test load pattern constant class.

@package test.model.load
@author ppizarror
@date 21/12/2018
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
#include "../../test_utils.h"
#include "../../../fnelem/model/loads/load_pattern_constant.h"

void __test_load_pattern_constant_init() {
    test_print_title("LOAD-PATTERN-CONSTANT", "load_pattern_constant");
    std::vector<Load *> *loads = new std::vector<Load *>();
    LoadPattern *ld = new LoadPatternConstant("LOAD", loads);
    ld->apply();
    ld->disp();
    delete ld;
    delete loads;
}

/**
 * Performs TEST-LOAD-PATTERN-CONSTANT suite.
 */
void test_load_pattern_constant_suite() {
    __test_load_pattern_constant_init();
}