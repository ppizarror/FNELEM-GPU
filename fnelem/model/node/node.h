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
#include "../../math/fematrix_utils.h"

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
    FEMatrix gdlid;

    // Coordinates of the node
    FEMatrix coords;

    // Displacements vector
    FEMatrix displ;

    // Loads vector
    FEMatrix loads;

    // Reaction vector
    FEMatrix reaction;

    // Init internal variables
    void init();

    // Destroy malloc inner variables
    void destroy();

    // Check external vector
    void check_vector(const FEMatrix *mat, std::string vector_name);

public:

    // Destroy node
    ~Node();

    // 2D node
    Node(std::string tag, double posx, double posy);

    // 3D node
    Node(std::string tag, double posx, double posy, double posz);

    // Assign
    Node &operator=(const Node &node);

    // Return number of degrees of freedom
    int get_ngdl() const;

    // Return node coordinates
    FEMatrix get_coordinates() const;

    // Get GDLID of node
    FEMatrix get_gdlid() const;

    // Get results of loads
    FEMatrix get_load_results() const;

    // Get node displacements
    FEMatrix get_displacements() const;

    // Get node reactions
    FEMatrix get_reactions() const;

    // Set node GDLID
    void set_gdlid(int local_id, int global_id);

    // Set node GDLID from vector/matrix
    void set_gdlid(const FEMatrix *gdl);

    // Set node displacements
    void set_displacement(int local_id, double d);

    // Set node displacements from vector/matrix
    void set_displacement(const FEMatrix *d);

    // Apply load to node
    void apply_load(const FEMatrix *load);

    // Adds element inner stress to reactions
    void apply_element_stress(const FEMatrix *sigma);

    // Display node information to console
    void disp() const;

};

#endif // MODEL_NODE