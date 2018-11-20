/**
FNELEM-GPU NODE ELEMENT - NODE DEFINITION
Structural nodes.

@package fnelem.model.base
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
#include "../base/ModelComponent.h"

/**
 * Node element.
 */
class Node : public ModelComponent {
private:

    // Number of degres of freedom
    int ngdl = 0;

    // Coordinates of the node
    double *coords;

    // ID of degrees of freedom
    int *gdlid;

    // Displacements vector
    double *displ;

    // Loads vector
    double *loads;

    // Reaction vector
    double *reaction;

    // Init internal variables
    void init();

public:

    // Empty node
    Node();

    // Destroy node
    ~Node();

    // 2D node
    Node(std::string tag, double posx, double posy);

    // 3D node
    Node(std::string tag, double posx, double posy, double posz);

};

/**
 * Creates an empty node.
 */
Node::Node() {
    this->coords = new double[0];
    this->init();
}

/**
 * Destroy node.
 */
Node::~Node() {
    delete[] coords;
    delete[] gdlid;
    delete[] displ;
    delete[] loads;
    delete[] reaction;
}

/**
 * Creates a 2D node.
 *
 * @param tag Tag of the node
 * @param posx Position x
 * @param posy Position y
 */
Node::Node(std::string tag, double posx, double posy) : ModelComponent(tag) {
    this->ngdl = 2;
    this->coords = new double[2];
    this->coords[0] = posx;
    this->coords[1] = posy;
    this->init();
}

/**
 * Creates a 2D node.
 *
 * @param tag Tag of the node
 * @param posx Position x
 * @param posy Position y
 */
Node::Node(std::string tag, double posx, double posy, double posz) : ModelComponent(tag) {
    this->ngdl = 3;
    this->coords = new double[3];
    this->coords[0] = posx;
    this->coords[1] = posy;
    this->coords[2] = posz;
    this->init();
}

/**
 * Init internal variables.
 */
void Node::init() {

    // Init ID of degrees of freedom
    this->gdlid = new int[this->ngdl];

    // Init displacements
    this->displ = new double[this->ngdl];

    // Init reactions
    this->reaction = new double[this->ngdl];

    // Init node loads
    this->loads = new double[this->ngdl];

    // Save initial values
    for (int i = 0; i < this->ngdl; i++) {
        this->gdlid[i] = -1;
        this->displ[i] = 0;
        this->reaction[i] = 0;
        this->loads[i] = 0;
    }

}
