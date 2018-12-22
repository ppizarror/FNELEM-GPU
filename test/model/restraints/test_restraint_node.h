/**
FNELEM-GPU - TEST
Test node restraints.

@package test.model.restraints
@author ppizarror
@date 30/11/2018
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
#include "../../../fnelem/model/restraints/restraint_node.h"

void __test_restraint_node() {
    test_print_title("RESTRAINT-NODE", "test_restraint_node");

    // Create node
    Node *n = new Node("NODE", 1.5, 3.2, 5.6);
    RestraintNode *r = new RestraintNode("R1", n);

    // Apply restraints, no DOFID has been setted, so GLOBALID of node is 0
    r->apply();
    FEMatrix *dofn1 = n->get_dof();
    assert(dofn1->is_double(0));
    r->disp();

    // Create restraints
    r->add_dofid(1);
    r->apply();
    assert(is_num_equal(n->get_dof(1), -1));
    assert(is_num_equal(n->get_dof(2), 0));
    assert(is_num_equal(n->get_dof(3), 0));

    // Apply more constraints
    r->add_dofid(3);
    r->apply();
    assert(is_num_equal(n->get_dof(1), -1));
    assert(is_num_equal(n->get_dof(2), 0));
    assert(is_num_equal(n->get_dof(3), -1));

    // Display restraint information
    r->disp();

    // Delete vars
    delete n;
    delete r;
    delete dofn1;

}

/**
 * Performs TEST-RESTRAINT-NODE suite.
 */
void test_restraint_node_suite() {
    __test_restraint_node();
}