/**
FNELEM-GPU - TEST LOAD
Test general load node.

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
#include "../../../fnelem/model/loads/load_node.h"

void __test_load_node_init() {

    // Creates node
    Node *n = new Node("NODE", 250, 0, -543);

    // Creates vector load
    FEMatrix *loadv = FEMatrix_vector(3);
    loadv->set(0, -5);

    // Creates load node
    LoadNode *ln = new LoadNode("LN", n, loadv);
    ln->disp();

    // Apply load node force
    n->disp();
    ln->apply(1.0);
    n->disp();

    // Var deletion
    delete n;
    delete loadv;
    delete ln;

}

/**
 * Performs TEST-LOAD-NODE suite.
 */
void test_load_node_suite() {
    __test_load_node_init();
}