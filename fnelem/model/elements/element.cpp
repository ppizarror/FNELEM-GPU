/**
FNELEM-GPU ELEMENTS
Base structural element class.

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
#include "element.h"

/**
 * Class destructor.
 */
Element::~Element() {
    delete this->nodes;
    if (this->initialized) {
        delete this->dofid;
        delete this->stiffness_global;
        delete this->stiffness_local;
        delete this->constitutive;
    }
}

/**
 * Return number of nodes.
 *
 * @return
 */
int Element::get_node_number() const {
    return this->nnodes;
}

/**
 * Default constructor.
 *
 * @param tag Element tag
 */
Element::Element(std::string tag) : ModelComponent(std::move(tag)) {}

/**
 * Get element nodes.
 *
 * @return
 */
std::vector<Node *> *Element::get_nodes() const {
    return this->nodes;
}

/**
 * Get number of degrees of freedom.
 *
 * @return
 */
int Element::get_ndof() const {
    return this->ndof;
}

/**
 * Get DOFID associated with the element, used by analysis method for create
 * matrix stiffness.
 *
 * @return
 */
FEMatrix *Element::get_dofid() const {
    return this->dofid->clone();
}

/**
 * Get local stiffness matrix.
 *
 * @return
 */
FEMatrix *Element::get_stiffness_local() const {
    return this->stiffness_local->clone();
}

/**
 * Get global stiffness matrix.
 *
 * @return
 */
FEMatrix *Element::get_stiffness_global() const {
    return this->stiffness_global->clone();
}

/**
 * Get local resistant force.
 *
 * @return
 */
FEMatrix *Element::get_force_local() const {
    return new FEMatrix(this->ndof, 1);
}

/**
 * Get global resistant force.
 *
 * @return
 */
FEMatrix *Element::get_force_global() const {
    return new FEMatrix(this->ndof, 1);
}

/**
 * Disp element information to console.
 */
void Element::disp() const {
    ModelComponent::disp();
}

/**
 * Returns constitutive matrix.
 *
 * @return
 */
FEMatrix *Element::get_constitutive() const {
    return this->constitutive->clone();
}

/**
 * Set as initialized.
 */
void Element::initialize() {
    this->initialized = true;
}