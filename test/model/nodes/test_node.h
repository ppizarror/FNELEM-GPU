/**
FNELEM-GPU - TEST
Test structure node class.

@package test.model.node
@author ppizarror
@date 19/11/2018
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
#include "../../../fnelem/model/nodes/node.h"

void __test_node_creation() {
    test_print_title("NODE", "test_node_creation");

    // Test 2D node
    Node *n = new Node("NODE1", 0, 1);
    n->disp();
    assert(n->get_model_tag() == "NODE1");
    FEMatrix *gdlid = new FEMatrix(2, 1);
    gdlid->set(0, 0);
    gdlid->set(1, 0);
    FEMatrix *id = n->get_dofid();
    assert(gdlid->equals(id));
    delete id;
    delete n;
    delete gdlid;

    // Test 3D node
    Node *n2 = new Node("NODE3D", 1.5, 3.2, 5.6);
    assert(n2->get_ndof() == 3);
    assert(is_num_equal(n2->get_pos_x(), 1.5));
    assert(is_num_equal(n2->get_pos_y(), 3.2));
    assert(is_num_equal(n2->get_pos_z(), 5.6));
    delete n2;

}

void __test_node_coordinates() {
    test_print_title("NODE", "test_coordinates");
    Node *n = new Node("NODE", 1, 2);
    FEMatrix *coords = n->get_coordinates();
    assert(coords->get(0) == 1);
    assert(coords->get(1) == 2);
    FEMatrix *id = n->get_dofid();
    assert(id->is_double(0));
    delete coords;
    delete id;
    delete n;
}

void __test_node_loads() {
    test_print_title("NODE", "test_loads");
    Node *n = new Node("NODE", 4, 5, -7);
    FEMatrix *loads = n->get_load_results();
    FEMatrix *displ = n->get_displacements();
    FEMatrix *react = n->get_reactions();
    assert(loads->is_zeros());
    assert(displ->is_zeros());
    assert(react->is_zeros());
    delete loads;
    delete displ;
    delete react;
    delete n;
}

void __test_node_set_gdlid() {
    test_print_title("NODE", "test_set_gdlid");
    Node *n1 = new Node("NODE", 0, 0, 0);
    FEMatrix *gdl = FEMatrix_vector(3);

    gdl->set(0, 5);
    gdl->set(1, 6);
    gdl->set(2, 8);
    n1->set_dof(1, 5);
    n1->set_dof(2, 6);
    n1->set_dof(3, 8);

    FEMatrix *id1 = n1->get_dofid();
    assert(*id1 == *gdl);
    delete n1;

    Node *n2 = new Node("NODE", 0, 0, 0);
    FEMatrix *id2 = n2->get_dofid();
    assert(id2->is_double(0));
    n2->set_dof(gdl);
    FEMatrix *id3 = n2->get_dofid();
    assert(id3->equals(gdl));
    delete n2;
    delete id1;
    delete id2;
    delete id3;
    delete gdl;
}

void __test_node_displacements() {
    test_print_title("NODE", "test_displacements");
    Node *n = new Node("NODE", 0, 0, 0);
    FEMatrix *displ = FEMatrix_vector(3);
    displ->set(0, 5);

    displ->set(1, -6);
    displ->set(2, 0);
    n->set_displacement(displ);
    FEMatrix *d1 = n->get_displacements();
    assert(*d1 == *displ);
    displ->set(0, -5);
    FEMatrix *d2 = n->get_displacements();
    assert(*d2 != *displ);
    n->set_displacement(1, -5);
    FEMatrix *d3 = n->get_displacements();
    assert(*d3 == *displ);

    delete d1;
    delete d2;
    delete d3;
    delete n;
    delete displ;
}

void __test_node_load() {
    test_print_title("NODE", "test_node_load");
    Node *n = new Node("NODE", 0, 0);
    FEMatrix *load = FEMatrix_vector(2);
    load->set(0, 5);
    load->set(1, -6);
    n->apply_load(load);
    FEMatrix *r1 = n->get_reactions();
    r1->disp();
    FEMatrix *minload = -*load;
    assert(r1->equals(minload));
    n->apply_element_stress(load);
    FEMatrix *r2 = n->get_reactions();
    assert(r2->is_zeros());
    delete r2;
    delete n;
    delete r1;
    delete load;
    delete minload;
}

void __test_node_full() {
    test_print_title("NODE", "test_node_full");
    Node *n = new Node("NODE", 1, 1, -1);
    FEMatrix *displ = FEMatrix_vector(3);
    displ->set(0, 5);
    displ->set(1, -6);
    displ->set(2, 0);
    n->set_displacement(displ);
    FEMatrix *load = FEMatrix_vector(3);
    load->set(0, 2);
    load->set(1, -1);
    n->apply_load(load);
    n->disp();
    delete n;
    delete displ;
    delete load;
}

void __test_node_simple() {
    test_print_title("NODE", "test_simple");
    Node *n1 = new Node("NODE", 1, 2);
    n1->disp();
    delete n1;
}

void __test_node_save_to_file() {
    test_print_title("NODE", "test_save_to_file");

    // Create node
    Node *n = new Node("N1", 1, 1, -1);
    FEMatrix *displ = FEMatrix_vector(3);
    displ->set(0, 5);
    displ->set(1, -6);
    displ->set(2, 0);
    n->set_displacement(displ);
    FEMatrix *load = FEMatrix_vector(3);
    load->set(0, 2);
    load->set(1, -1);
    n->apply_load(load);

    // Create file
    std::ofstream plik;
    plik.open("out/test-node-properties.txt");
    n->save_properties(plik);
    n->save_displacements(plik);
    n->save_reactions(plik);
    plik.close();

    delete n;
    delete displ;
    delete load;
}

/**
 * Performs TEST-NODE suite.
 */
void test_node_suite() {
    __test_node_simple();
    __test_node_creation();
    __test_node_coordinates();
    __test_node_loads();
    __test_node_set_gdlid();
    __test_node_displacements();
    __test_node_load();
    __test_node_full();
    __test_node_save_to_file();
}