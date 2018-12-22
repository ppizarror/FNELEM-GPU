/**
FNELEM-GPU - TEST MEMBRANE
Test membrane element.

@package test.model.elements
@author ppizarror
@date 26/11/2018
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
#include "../../../fnelem/model/elements/membrane.h"

void __test_membrane_creation() {
    test_print_title("ELEMENTS-MEMBRANE", "test_membrane_creation");

    // Node creation
    Node *n1 = new Node("N1", 0, 0);
    Node *n2 = new Node("N2", 250, 0);
    Node *n3 = new Node("N3", 250, 100);
    Node *n4 = new Node("N4", 0, 100);

    // Creates membrane
    Membrane *mem = new Membrane("MEM", n1, n2, n3, n4, 300000, 0.15, 20);

    // Check forces are zero
    FEMatrix *localf = mem->get_force_local();
    FEMatrix *globalf = mem->get_force_global();

    assert(localf->is_zeros());
    assert(globalf->is_zeros());

    // Check node pointers are the same
    std::vector<Node *> *nodes = mem->get_nodes();
    assert (n1 == nodes->at(0));
    assert (n2 == nodes->at(1));
    assert (n3 == nodes->at(2));
    assert (n4 == nodes->at(3));

    // Test dimension
    assert(is_num_equal(mem->get_width(), 250));
    assert(is_num_equal(mem->get_height(), 100));

    // Test ndof
    assert(is_num_equal(mem->get_node_number(), 4));
    assert(is_num_equal(mem->get_ndof(), 8));

    // Display membrane information
    mem->disp();

    // Test constitutive matrix
    FEMatrix *c = mem->get_constitutive();
    (*c) *= pow(10, -5);
    assert(is_num_equal(c->get(0, 0), 3.0690537084398981094));
    assert(is_num_equal(c->get(0, 1), 0.46035805626598458318));
    assert(is_num_equal(c->get(0, 2), 0));
    assert(is_num_equal(c->get(1, 0), 0.46035805626598458318));
    assert(is_num_equal(c->get(1, 1), 3.0690537084398981094));
    assert(is_num_equal(c->get(1, 2), 0));
    assert(is_num_equal(c->get(2, 0), 0));
    assert(is_num_equal(c->get(2, 1), 0));
    assert(is_num_equal(c->get(2, 2), 1.3043478260869567631));

    // Test local stiffness of matrix
    FEMatrix *k = mem->get_stiffness_local();
    FEMatrix *kl = mem->get_stiffness_local();
    FEMatrix *ke = mem->get_stiffness_global();
    (*k) *= pow(10, -6);

    assert(k->is_symmetric());
    assert(k->is_square());

    assert(is_num_equal(k->get(0, 0), 2.9923273657289004568));
    assert(is_num_equal(k->get(0, 1), 0.88235294117647056211));
    assert(is_num_equal(k->get(0, 2), 0.26854219948849117339));
    assert(is_num_equal(k->get(0, 3), -0.42199488491048597893));
    assert(is_num_equal(k->get(0, 4), -1.4961636828644502284));
    assert(is_num_equal(k->get(0, 5), -0.88235294117647056211));
    assert(is_num_equal(k->get(0, 6), -1.7647058823529413463));
    assert(is_num_equal(k->get(0, 7), 0.42199488491048597893));
    assert(is_num_equal(k->get(1, 1), 7.2890025575447561224));
    assert(is_num_equal(k->get(1, 2), 0.42199488491048597893));
    assert(is_num_equal(k->get(1, 3), 2.2097186700767257328));
    assert(is_num_equal(k->get(1, 4), -0.88235294117647056211));
    assert(is_num_equal(k->get(1, 5), -2.731457800511508438));
    assert(is_num_equal(k->get(1, 6), -0.42199488491048597893));
    assert(is_num_equal(k->get(1, 7), -4.9411764705882346149));
    assert(is_num_equal(k->get(2, 2), 2.9923273657289004568));
    assert(is_num_equal(k->get(2, 3), -0.88235294117647056211));
    assert(is_num_equal(k->get(2, 4), -1.7647058823529413463));
    assert(is_num_equal(k->get(2, 5), -0.42199488491048597893));
    assert(is_num_equal(k->get(2, 6), -1.4961636828644502284));
    assert(is_num_equal(k->get(2, 7), 0.88235294117647056211));
    assert(is_num_equal(k->get(3, 3), 7.2890025575447561224));
    assert(is_num_equal(k->get(3, 4), 0.42199488491048597893));
    assert(is_num_equal(k->get(3, 5), -4.9411764705882346149));
    assert(is_num_equal(k->get(3, 6), 0.88235294117647056211));
    assert(is_num_equal(k->get(3, 7), -2.731457800511508438));
    assert(is_num_equal(k->get(4, 4), 2.9923273657289004568));
    assert(is_num_equal(k->get(4, 5), 0.88235294117647056211));
    assert(is_num_equal(k->get(4, 6), 0.26854219948849117339));
    assert(is_num_equal(k->get(4, 7), -0.42199488491048597893));
    assert(is_num_equal(k->get(5, 5), 7.2890025575447561224));
    assert(is_num_equal(k->get(5, 6), 0.42199488491048597893));
    assert(is_num_equal(k->get(5, 7), 2.2097186700767257328));
    assert(is_num_equal(k->get(6, 6), 2.9923273657289004568));
    assert(is_num_equal(k->get(6, 7), -0.88235294117647056211));
    assert(is_num_equal(k->get(7, 7), 7.2890025575447561224));

    // Check local and global matrices are the same
    assert(kl->equals(ke));

    // Test local degrees of freedom, must be 0
    mem->set_dofid();
    FEMatrix *dofid = mem->get_dofid();
    assert(dofid->is_double(0));

    // Delete vars
    delete n1;
    delete n2;
    delete n3;
    delete n4;
    delete mem;
    delete k;
    delete kl;
    delete ke;
    delete localf;
    delete globalf;
    delete dofid;
    delete c;

}

void __test_membrane_evalxy() {
    test_print_title("ELEMENTS-MEMBRANE", "test_membrane_evalxy");

    // Node creation
    Node *n1 = new Node("N1", 0, 0);
    Node *n2 = new Node("N2", 250, 0);
    Node *n3 = new Node("N3", 250, 100);
    Node *n4 = new Node("N4", 0, 100);

    // Creates membrane
    Membrane *mem = new Membrane("MEM", n1, n2, n3, n4, 300000, 0.15, 20);

    // Calcule displacement
    FEMatrix *d = mem->get_displacement(0, 0);
    assert(d->length() == 2);
    assert(d->is_zeros());

    // Calculate deformation/strain
    FEMatrix *st = mem->get_deformation(0, 0);
    assert(st->length() == 3);
    assert(st->is_zeros());

    // Calculate stress
    FEMatrix *sr = mem->get_stress(0, 0);
    assert(sr->length() == 3);
    assert(sr->is_zeros());

    // Var deletion
    delete n1;
    delete n2;
    delete n3;
    delete n4;
    delete mem;
    delete d;
    delete st;
    delete sr;

}

void __test_membrane_forces() {
    test_print_title("ELEMENTS-MEMBRANE", "test_membrane_forces");

    // Node creation
    Node *n1 = new Node("N1", 0, 0);
    Node *n2 = new Node("N2", 250, 0);
    Node *n3 = new Node("N3", 250, 100);
    Node *n4 = new Node("N4", 0, 100);

    // Creates membrane
    Membrane *mem = new Membrane("MEM", n1, n2, n3, n4, 300000, 0.15, 20);

    // Add force
    FEMatrix *f = FEMatrix_vector(2);
    f->set(0, -1);
    f->set(1, 5);
    mem->add_equivalent_force_node(4, f);
    mem->add_equivalent_force_node(4, f);

    // Get local/global force
    FEMatrix *fl = mem->get_force_local();
    FEMatrix *fg = mem->get_force_global();
    assert(fl->is_zeros());
    assert(is_num_equal(fg->get(6), 2));
    assert(is_num_equal(fg->get(7), -10)); // Fglobal = Flocal - Feq

    // Add force to reactions
    mem->add_force_to_reaction();
    assert(is_num_equal(n1->get_reaction(1), 0));
    assert(is_num_equal(n1->get_reaction(2), 0));
    assert(is_num_equal(n2->get_reaction(1), 0));
    assert(is_num_equal(n2->get_reaction(2), 0));
    assert(is_num_equal(n3->get_reaction(1), 0));
    assert(is_num_equal(n3->get_reaction(2), 0));
    assert(is_num_equal(n4->get_reaction(1), 0));
    assert(is_num_equal(n4->get_reaction(2), 0));

    // Variable deletion
    delete n1;
    delete n2;
    delete n3;
    delete n4;
    delete mem;
    delete f;
    delete fl;
    delete fg;

}

void __test_membrane_save_to_file() {
    test_print_title("ELEMENTS-MEMBRANE", "test_membrane_save_to_file");

    // Node creation
    Node *n1 = new Node("N1", 0, 0);
    Node *n2 = new Node("N2", 250, 0);
    Node *n3 = new Node("N3", 250, 100);
    Node *n4 = new Node("N4", 0, 100);

    // Creates membrane
    Membrane *mem = new Membrane("MEM", n1, n2, n3, n4, 300000, 0.15, 20);

    // Save file
    std::ofstream plik;
    plik.open("out/test-save-membrane.txt");
    mem->save_properties(plik);
    plik << "\n" << std::endl;
    mem->save_internal_stress(plik);
    plik.close();

    // Delete vars
    delete n1;
    delete n2;
    delete n3;
    delete n4;
    delete mem;

}

/**
 * Performs TEST-MEMBRANE suite.
 */
void test_membrane_suite() {
    __test_membrane_creation();
    __test_membrane_evalxy();
    __test_membrane_forces();
    __test_membrane_save_to_file();
}