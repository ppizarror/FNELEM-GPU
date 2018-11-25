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
}

void test_loads() {
    Node n = Node("NODE", 4, 5, -7);
    assert(n.get_load_results().is)
}

int main() {
    test_node_creation();
    test_coordinates();
    test_loads();
    return 0;
}