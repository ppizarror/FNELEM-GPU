/**
FNELEM-GPU NODE ELEMENT - NODE DEFINITION
Structural nodes.

@package fnelem.model.node
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

// Include header
#include "node.h"

/**
 * Destroy inner variables.
 */
void Node::destroy() {
}

/**
 * Destroy node.
 */
Node::~Node() {
    this->destroy();
};

/**
 * Creates a 2D node.
 *
 * @param tag Tag of the node
 * @param posx Position x
 * @param posy Position y
 */
Node::Node(std::string tag, double posx, double posy) : ModelComponent(std::move(tag)) {
    this->ngdl = 2;
    this->coords = FEMatrix_vector(this->ngdl);
    this->coords.set(0, posx);
    this->coords.set(1, posy);
    this->init();
}

/**
 * Creates a 2D node.
 *
 * @param tag Tag of the node
 * @param posx Position x
 * @param posy Position y
 */
Node::Node(std::string tag, double posx, double posy, double posz) : ModelComponent(std::move(tag)) {
    this->ngdl = 3;
    this->coords = FEMatrix_vector(this->ngdl);
    this->coords.set(0, posx);
    this->coords.set(1, posy);
    this->coords.set(2, posz);
    this->init();
}

/**
 * Init internal variables.
 */
void Node::init() {

    // Init ID of degrees of freedom
    this->gdlid = FEMatrix_vector(this->ngdl);

    // Init displacements
    this->displ = FEMatrix_vector(this->ngdl);

    // Init reactions
    this->reaction = FEMatrix_vector(this->ngdl);

    // Init node loads
    this->loads = FEMatrix_vector(this->ngdl);

    // Save initial values
    for (int i = 0; i < this->ngdl; i++) {
        this->gdlid.set(i, -1);
        this->displ.set(i, 0);
        this->reaction.set(i, 0);
        this->loads.set(i, 0);
    }

}

/**
 * Return number of degrees of freedom.
 *
 * @return Ngdl
 */
int Node::get_ngdl() const {
    return this->ngdl;
}

/**
 * Return a copy of node coordinates.
 *
 * @return Coordinates
 */
FEMatrix Node::get_coordinates() const {
    return this->coords.clone();
}

/**
 * Get GDLID vector.
 *
 * @return
 */
FEMatrix Node::get_gdlid() const {
    return this->gdlid.clone();
}

/**
 * Assign operator.
 *
 * @param node
 * @return
 */
Node &Node::operator=(const Node &node) {

    // Destroy inner vectors
    this->destroy();

    // Assign operator
    this->ngdl = node.ngdl;
    this->coords = node.coords;
    this->displ = node.displ;
    this->reaction = node.reaction;
    this->loads = node.loads;

    // Call ModelComponent assign
    ModelComponent::operator=(node);
    return *this;

}