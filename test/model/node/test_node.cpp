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
#include "../../../fnelem/model/node/node.h"

void test_node_creation() {

    // Test nodo en 2D
    Node n = Node("NODE1", 0, 1);
    assert(n.get_model_tag() == "NODE1");
    FEMatrix gdlid = FEMatrix(2, 1);
    gdlid.set(0, -1);
    gdlid.set(1, -1);
    assert(gdlid == n.get_gdlid());

    // Test nodo en 3D
    n = Node("NODE3D", 1.5, 3.2, 5.6);
    assert(n.get_ngdl() == 3);

}

void test_coordinates() {
    Node n = Node("NODE", 1, 2);
    assert(n.get_coordinates().get(0) == 1);
    assert(n.get_coordinates().get(1) == 2);
    assert(n.get_gdlid().is_double(-1));
}

void test_loads() {
    Node n = Node("NODE", 4, 5, -7);
    assert(n.get_load_results().is_zeros());
    assert(n.get_displacements().is_zeros());
    assert(n.get_reactions().is_zeros());
}

void test_set_gdlid() {
    Node n = Node("NODE", 0, 0, 0);
    FEMatrix gdl = FEMatrix_vector(3);
    gdl.set(0, 5);
    gdl.set(1, 6);
    gdl.set(2, 8);
    n.set_gdlid(1, 5);
    n.set_gdlid(2, 6);
    n.set_gdlid(3, 8);
    assert(n.get_gdlid() == gdl);
    n = Node("NODE", 0, 0, 0);
    assert(n.get_gdlid().is_double(-1));
    n.set_gdlid(&gdl);
    assert(n.get_gdlid() == gdl);
}

void test_node_displacements() {
    Node n = Node("NODE", 0, 0, 0);
    FEMatrix displ = FEMatrix_vector(3);
    displ.set(0, 5);
    displ.set(1, -6);
    displ.set(2, 0);
    n.set_displacement(&displ);
    assert(n.get_displacements() == displ);
    displ.set(0, -5);
    assert(n.get_displacements() != displ);
    n.set_displacement(1, -5);
    assert(n.get_displacements() == displ);
}

void test_node_load() {
    Node n = Node("NODE", 0, 0);
    FEMatrix load = FEMatrix_vector(2);
    load.set(0, 5);
    load.set(1, -6);
    n.apply_load(&load);
    assert(n.get_reactions() == -load);
    n.apply_element_stress(&load);
    assert(n.get_reactions().is_zeros());
}

void test_node_full() {
    Node n = Node("NODE", 0, 0, 0);
    FEMatrix displ = FEMatrix_vector(3);
    displ.set(0, 5);
    displ.set(1, -6);
    displ.set(2, 0);
    n.set_displacement(&displ);
    FEMatrix load = FEMatrix_vector(3);
    load.set(0, 2);
    load.set(1, -1);
    n.apply_load(&load);
    n.disp();
}

int main() {
    test_node_creation();
    test_coordinates();
    test_loads();
    test_set_gdlid();
    test_node_displacements();
    test_node_load();
    test_node_full();
    return 0;
}