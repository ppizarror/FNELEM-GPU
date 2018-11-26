#include <utility>

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

// Include source
#include "element.h"

/**
 * Class destructor.
 */
Element::~Element() {
    delete this->nodes;
    if (this->initialized) {
        delete this->gdlid;
        delete this->force_global;
        delete this->force_local;
        delete this->stiffness_global;
        delete this->stiffness_local;
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
int Element::get_ngdl() const {
    return this->ngdl;
}

/**
 * Get GDLID associated with the element.
 *
 * @return
 */
FEMatrix *Element::get_gdlid() const {
    return this->gdlid->clone();
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
    return this->force_local->clone();
}

/**
 * Get global resistant force.
 *
 * @return
 */
FEMatrix *Element::get_force_global() const {
    return this->force_global->clone();
}

/**
 * Disp element information to console.
 */
void Element::disp() const {
    ModelComponent::disp();
}