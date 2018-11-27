#include <utility>

/**
FNELEM-GPU MEMBRANE ELEMENT
Bidimensional membrane element composed by 4 nodes.

@package fnelem.model.elements
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

// Include header
#include "membrane.h"

/**
 * Membrane constructor.
 *
 * @param tag Membrane tag
 * @param n1 Node 1
 * @param n2 Node 2
 * @param n3 Node 3
 * @param n4 Node 4
 * @param E Elasticity modulus of the section
 * @param poisson Poisson constant of the section
 * @param thickness Thickness of the section
 */
Membrane::Membrane(std::string tag, Node *n1, Node *n2, Node *n3, Node *n4, double E, double poisson, double thickness)
        : Element(std::move(tag)) {

    // Stores material
    this->ngdl = 8;
    this->E = E;
    this->poisson = poisson;

    // Stores nodes
    this->nodes->push_back(n1);
    this->nodes->push_back(n2);
    this->nodes->push_back(n3);
    this->nodes->push_back(n4);

    // Generate constitutive matrix
    this->constitutive = new FEMatrix(3, 3);
    this->constitutive->set(0, 0, 1 / (1 - poisson * poisson));
    this->constitutive->set(0, 1, poisson / (1 - poisson * poisson));
    this->constitutive->set(1, 0, poisson / (1 - poisson * poisson));
    this->constitutive->set(1, 1, 1 / (1 - poisson * poisson));
    this->constitutive->set(2, 2, 1 / (2 + 2 * poisson));

    // Calculates dimension
    //     4 ------------- 3
    //     |       y       |
    //  2h |       #x      |
    //     |               |
    //     1 ------------- 2
    //            2b
    double db1, db2, dh1, dh2;
    db1 = fabs(n1->get_pos_x() - n2->get_pos_x()) / 2;
    db2 = fabs(n3->get_pos_x() - n4->get_pos_x()) / 2;
    dh1 = fabs(n1->get_pos_y() - n4->get_pos_y()) / 2;
    dh2 = fabs(n2->get_pos_y() - n3->get_pos_y()) / 2;

    if (fabs(db1 - db2) > __FEMATRIX_ZERO_TOL) {
        throw std::logic_error("[MEMBRANE] Invalid node dimension at @x");
    }
    if (fabs(dh1 - dh2) > __FEMATRIX_ZERO_TOL) {
        throw std::logic_error("[MEMBRANE] Invalid node dimension at @y");
    }

    // Stores geometry
    this->t = thickness;
    this->b = db1;
    this->h = dh1;

    // Init forces
    this->Feq = FEMatrix_vector(8);

    // Init matrices
    this->Feq = new FEMatrix(8, 8);

}

/**
 * Destructor.
 */
Membrane::~Membrane() {
    delete this->Feq;
}