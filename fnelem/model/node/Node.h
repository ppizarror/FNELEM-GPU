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

// Init header file
#ifndef MODEL_NODE_HEADER_FILE
#define MODEL_NODE_HEADER_FILE

// Include headers
#include "../base/model_component.h"
#include "../../math/fematrix.h"

// Library imports
#include <iostream>

/**
 * Node element.
 */
class Node : public ModelComponent {
private:

    // Number of degres of freedom
    int ngdl = 0;

    // ID of degrees of freedom
    int *gdlid;

    // Coordinates of the node
    double *coords;

    // Displacements vector
    double *displ;

    // Loads vector
    double *loads;

    // Reaction vector
    double *reaction;

    // Init internal variables
    void init();

public:

    // Destroy node
    ~Node();

    // 2D node
    Node(std::string tag, double posx, double posy);

    // 3D node
    Node(std::string tag, double posx, double posy, double posz);

    // Return number of degrees of freedom
    int get_ngdl() const;

    // Return node coordinates
    double *get_coordinates() const;

};

#endif // MODEL_NODE