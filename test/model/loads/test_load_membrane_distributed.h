/**
FNELEM-GPU - TEST LOAD
Test membrane distributed load.

@package test.model.load
@author ppizarror
@date 01/12/2018
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
#include "../../../fnelem/model/loads/load_membrane_distributed.h"

void __test_load_membrane_distributed_init() {
    test_print_title("LOAD-MEMBRANE-DISTRIBUTED", "test_load_membrane_distributed_init");

    // Creates membrane
    Node *n1 = new Node("N1", 0, 0);
    Node *n2 = new Node("N2", 250, 0);
    Node *n3 = new Node("N3", 250, 100);
    Node *n4 = new Node("N4", 0, 100);
    Membrane *mem = new Membrane("MEM", n1, n2, n3, n4, 300000, 0.15, 20);

    // Creates load
    LoadMembraneDistributed *load = new LoadMembraneDistributed("LOAD-MEM", mem, 1, 2, 10, 0, 10, 1);
    load->apply(1.0);
    load->disp();

    // Var deletion
    delete mem;
    delete n1;
    delete n2;
    delete n3;
    delete n4;
    delete load;

}

/**
 * Performs TEST-MEMBRANED-DISTRIBUTED suite.
 */
void test_load_membrane_distributed_suite() {
    __test_load_membrane_distributed_init();
}